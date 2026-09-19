from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")
pytest.importorskip("pytestqt")

from PySide6.QtCore import QPointF

from qss_telemetry.fixtures import make_demo_result
from qss_telemetry.schema import AxisData, ChannelData
from qss_telemetry.widgets.plot_stack import PlotStack


def _lap_with_channel_values(lap, channel_id, values):
    channel = lap.channels[channel_id]
    values = np.asarray(values, dtype=float)
    updated = ChannelData(channel.metadata, values, np.ones(values.shape, dtype=bool))
    return replace(lap, channels={**lap.channels, channel_id: updated})


def _single_channel_lap(template, axis, values):
    channel = template.channels["speed_mps"]
    values = np.asarray(values, dtype=float)
    updated = ChannelData(channel.metadata, values, np.ones(values.shape, dtype=bool))
    return replace(template, axes={"native": axis}, channels={"speed_mps": updated})


def test_dragging_any_cursor_line_moves_the_shared_cursor(qtbot):
    stack = PlotStack()
    qtbot.addWidget(stack)
    lap = make_demo_result().get_lap("baseline", "flying")
    stack.set_data(lap, ["speed_mps", "gear"], cursor_value=1.0)
    values = []
    stack.cursorValueChanged.connect(values.append)

    stack.cursor_lines[1].setValue(2.25)

    assert all(line.value() == pytest.approx(2.25) for line in stack.cursor_lines)
    assert values[-1] == pytest.approx(2.25)


def test_waveforms_disable_mouse_pan_and_keep_linked_toolbar_zoom(qtbot):
    stack = PlotStack()
    qtbot.addWidget(stack)
    lap = make_demo_result().get_lap("baseline", "flying")
    stack.set_data(lap, ["speed_mps", "gear"], cursor_value=0.0)

    for plot in stack.plots:
        assert list(plot.getViewBox().state["mouseEnabled"]) == [True, False]
        assert plot.getViewBox().pan_enabled is False

    stack.fit_x_range()
    initial = stack.plots[0].getViewBox().viewRange()[0]
    stack.zoom_in()
    zoomed = stack.plots[0].getViewBox().viewRange()[0]
    assert (zoomed[1] - zoomed[0]) < (initial[1] - initial[0])
    assert all(plot.getViewBox().viewRange()[0] == pytest.approx(zoomed) for plot in stack.plots)

    stack.zoom_out()
    restored = stack.plots[0].getViewBox().viewRange()[0]
    assert (restored[1] - restored[0]) > (zoomed[1] - zoomed[0])


def test_wheel_zoom_is_preserved_when_drag_pan_is_blocked(qtbot):
    stack = PlotStack()
    qtbot.addWidget(stack)
    lap = make_demo_result().get_lap("baseline", "flying")
    stack.set_data(lap, ["speed_mps", "gear"], cursor_value=0.0)
    view_box = stack.plots[0].getViewBox()
    before = view_box.viewRange()[0]

    class Wheel:
        def delta(self):
            return 120

        def pos(self):
            return QPointF(100.0, 100.0)

        def accept(self):
            pass

        def ignore(self):
            raise AssertionError("wheel zoom was disabled")

    view_box.wheelEvent(Wheel())

    after = view_box.viewRange()[0]
    assert (after[1] - after[0]) < (before[1] - before[0])


def test_delta_mode_draws_active_minus_datum_on_the_active_axis(qtbot):
    stack = PlotStack()
    qtbot.addWidget(stack)
    datum = make_demo_result().get_lap("baseline", "flying")
    active = make_demo_result().get_lap("baseline", "flying")
    active_speed = active.channels["speed_mps"]
    from dataclasses import replace

    from qss_telemetry.schema import ChannelData

    active = replace(
        active,
        channels={
            **active.channels,
            "speed_mps": ChannelData(active_speed.metadata, active_speed.values + 2.0, active_speed.valid),
        },
    )

    stack.set_data(active, ["speed_mps"], datum_lap=datum, display_mode="delta")

    assert stack.curves[("speed_mps", "delta")].getData()[1].tolist() == pytest.approx([2.0, 2.0, 2.0, 2.0])
    assert stack.curves.get(("speed_mps", "reference")) is None
    assert stack.curves.get(("speed_mps", "comparison")) is None


def test_delta_mode_renders_active_minus_datum_with_clear_label(qtbot):
    datum = make_demo_result().get_lap("baseline", "flying")
    active = _lap_with_channel_values(datum, "speed_mps", [12.0, 9.0, 14.0, 11.0])
    stack = PlotStack()
    qtbot.addWidget(stack)

    stack.set_data(active, ["speed_mps"], datum_lap=datum, display_mode="delta")

    assert stack.plot_count == 1
    assert set(stack.curves) == {("speed_mps", "delta")}
    _, values = stack.curves[("speed_mps", "delta")].getData()
    np.testing.assert_allclose(values, [2.0, -1.0, 4.0, 1.0])
    title = stack.plots[0].getPlotItem().titleLabel.text
    assert "Delta" in title
    assert "active - datum" in title
    assert stack.curves[("speed_mps", "delta")].opts["name"] == "Delta (active - datum)"


def test_delta_mode_aligns_active_and_datum_before_subtracting(qtbot):
    template = make_demo_result().get_lap("baseline", "flying")
    active = _single_channel_lap(
        template,
        AxisData("native", np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.0, 10.0, 20.0, 30.0])),
        [12.0, 12.0, 12.0, 12.0],
    )
    datum = _single_channel_lap(
        template,
        AxisData("native", np.array([0.0, 1.0, 2.0]), np.array([0.0, 15.0, 30.0])),
        [10.0, 20.0, 30.0],
    )
    stack = PlotStack()
    qtbot.addWidget(stack)

    stack.set_data(active, ["speed_mps"], datum_lap=datum, display_mode="delta", axis_mode="distance")

    x_values, delta_values = stack.curves[("speed_mps", "delta")].getData()
    np.testing.assert_allclose(x_values, [0.0, 10.0, 15.0, 20.0, 30.0])
    np.testing.assert_allclose(delta_values, [2.0, -14.0 / 3.0, -8.0, -34.0 / 3.0, -18.0])


def test_delta_mode_skips_channels_missing_from_datum_without_blocking_others(qtbot):
    datum = make_demo_result().get_lap("baseline", "flying")
    active = _lap_with_channel_values(datum, "speed_mps", [12.0, 12.0, 12.0, 12.0])
    datum_without_speed = replace(datum, channels={channel_id: channel for channel_id, channel in datum.channels.items() if channel_id != "speed_mps"})
    stack = PlotStack()
    qtbot.addWidget(stack)

    stack.set_data(active, ["speed_mps", "gear"], datum_lap=datum_without_speed, display_mode="delta")

    assert stack.plot_count == 1
    assert set(stack.curves) == {("gear", "delta")}


def test_delta_mode_keeps_shared_cursor_zoom_and_pan_behavior(qtbot):
    datum = make_demo_result().get_lap("baseline", "flying")
    active = _lap_with_channel_values(datum, "speed_mps", [12.0, 12.0, 12.0, 12.0])
    stack = PlotStack()
    qtbot.addWidget(stack)
    stack.set_data(active, ["speed_mps", "gear"], datum_lap=datum, display_mode="delta", cursor_value=1.0)

    assert all(list(plot.getViewBox().state["mouseEnabled"]) == [True, False] for plot in stack.plots)
    stack.cursor_lines[1].setValue(2.25)
    assert all(line.value() == pytest.approx(2.25) for line in stack.cursor_lines)

    initial = stack.plots[0].getViewBox().viewRange()[0]
    stack.zoom_in()
    zoomed = stack.plots[0].getViewBox().viewRange()[0]
    assert (zoomed[1] - zoomed[0]) < (initial[1] - initial[0])
    assert all(plot.getViewBox().viewRange()[0] == pytest.approx(zoomed) for plot in stack.plots)
