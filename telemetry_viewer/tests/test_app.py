"""Behavioral tests for the injected-result desktop viewer.

These tests intentionally use the real Qt and pyqtgraph widgets.  They are
skipped as a group when the optional desktop test dependencies are absent from
the active Python environment.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")
pytest.importorskip("pytestqt")

from qss_telemetry.app import TelemetryViewerWindow, ViewerModel
from qss_telemetry.fixtures import make_demo_result
from qss_telemetry.schema import AxisData, CaseData, ChannelData, LapData, ResultFile


def _result_with_comparison_lap() -> ResultFile:
    result = make_demo_result()
    reference = result.get_lap("baseline", "flying")
    slow_axis = AxisData("native", np.array([0.0, 1.25, 2.5, 3.75]), reference.axes["native"].distance_m)
    slow_channels = {
        channel_id: ChannelData(channel.metadata, channel.values.copy(), channel.valid.copy())
        for channel_id, channel in reference.channels.items()
    }
    comparison = LapData("slow", {"role": "comparison"}, {"native": slow_axis}, slow_channels, reference.track_id)
    case = replace(result.cases["baseline"], laps={"flying": reference, "slow": comparison})
    return replace(result, cases={"baseline": case})


def _result_with_speed_values(result: ResultFile, values: list[float]) -> ResultFile:
    reference = result.get_lap("baseline", "flying")
    speed = reference.channels["speed_mps"]
    updated_speed = ChannelData(speed.metadata, np.asarray(values, dtype=float), speed.valid.copy())
    lap = replace(reference, channels={**reference.channels, "speed_mps": updated_speed})
    case = replace(result.cases["baseline"], laps={"flying": lap})
    return replace(result, cases={"baseline": case})


def test_fixture_result_is_visible_in_browser_and_default_channels(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)

    assert window.case_lap_tree.topLevelItemCount() == 1
    assert window.case_lap_tree.topLevelItem(0).text(0) == "Baseline"
    assert window.channel_selector.is_checked("speed_mps")
    assert window.channel_selector.is_checked("long_accel_mps2")
    assert window.channel_selector.is_checked("lat_accel_mps2")
    assert window.channel_selector.is_checked("gear")


def test_model_can_select_a_loaded_result_as_datum_and_enable_delta_mode():
    datum = make_demo_result()
    active = _result_with_speed_values(datum, [12.0, 12.0, 12.0, 12.0])
    model = ViewerModel({"datum": datum, "active": active})

    model.select_result("active")
    model.set_datum_result("datum")
    model.set_display_mode("delta")

    assert model.datum_result_name == "datum"
    assert model.datum_lap is not None
    assert model.datum_lap.id == "flying"
    assert model.display_mode == "delta"


def test_window_renders_selected_result_delta_to_datum(qtbot):
    datum = make_demo_result()
    active = _result_with_speed_values(datum, [12.0, 12.0, 12.0, 12.0])
    window = TelemetryViewerWindow({"datum": datum, "active": active})
    qtbot.addWidget(window)

    window.select_result("active")
    window.set_datum_result("datum")
    window.set_display_mode("delta")

    curve = window.plot_stack.curves[("speed_mps", "delta")]
    assert curve.getData()[1].tolist() == pytest.approx([2.0, 2.0, 2.0, 2.0])
    assert window.waveform_mode_combo.currentData() == "delta"


def test_catalog_channels_are_included_in_default_view_when_available():
    result = make_demo_result()
    reference = result.get_lap("baseline", "flying")
    extra = {}
    for channel_id in ("aSteer", "SCL", "SCD", "throttle", "brake", "front_ride_height_in", "Fz_1"):
        metadata = replace(reference.channels["gear"].metadata, id=channel_id, label=channel_id)
        extra[channel_id] = ChannelData(metadata, np.zeros(4), np.ones(4, dtype=bool))
    lap = replace(reference, channels={**reference.channels, **extra})
    case = replace(result.cases["baseline"], laps={"flying": lap})

    model = ViewerModel(replace(result, cases={"baseline": case}))

    assert {"aSteer", "SCL", "SCD", "throttle", "brake", "front_ride_height_in", "Fz_1"} <= set(model.selected_channel_ids)
    assert "Aero" in model.grouped_channels()
    assert "Tire constraints" in model.grouped_channels()


def test_cursor_updates_all_plots_and_schematic_map_marker(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)

    window.set_cursor_index(2)

    assert window.cursor_index == 2
    assert window.map_widget.marker_distance_m == pytest.approx(20.0)
    assert window.map_widget.quality_label.text() == "Schematic map"
    assert all(line.value() == pytest.approx(2.0) for line in window.plot_stack.cursor_lines)


def test_dragging_waveform_cursor_updates_model_slider_and_map(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)

    window.plot_stack.cursor_lines[0].setValue(2.1)

    assert window.cursor_index == 2
    assert window.cursor_slider.value() == 2
    assert window.map_widget.marker_distance_m == pytest.approx(20.0)
    assert all(line.value() == pytest.approx(2.0) for line in window.plot_stack.cursor_lines)


def test_fullscreen_toggle_restores_the_previous_window_size(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)
    window.resize(1111, 777)
    window.show()
    qtbot.wait(20)
    normal_size = window.size()

    window.set_fullscreen(True)
    qtbot.wait(20)
    assert window.isFullScreen()

    window.set_fullscreen(False)
    qtbot.wait(20)
    assert not window.isFullScreen()
    assert window.size() == normal_size


def test_initial_window_fits_available_screen(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)
    window.show()
    qtbot.wait(20)

    available = window.screen().availableGeometry()
    assert window.width() <= available.width()
    assert window.height() <= available.height()


def test_fullscreen_stays_fullscreen_when_slider_moves_slap(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)
    window.show()
    qtbot.wait(20)

    window.set_fullscreen(True)
    qtbot.wait(20)
    assert window.isFullScreen()

    window.cursor_slider.setValue(1)
    qtbot.wait(20)

    assert window.isFullScreen()


def test_manual_slap_slider_advances_cursor_and_map_deterministically(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)
    window.set_cursor_index(0)
    window.cursor_slider.setValue(1)

    assert window.cursor_index == 1
    assert window.map_widget.marker_distance_m == pytest.approx(10.0)


def test_viewer_keeps_manual_slap_slider_without_playback_controls(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)

    assert not hasattr(window, "play_button")
    assert not hasattr(window, "speed_combo")
    assert not hasattr(window, "playback_controller")

    window.cursor_slider.setValue(2)

    assert window.cursor_index == 2


def test_absent_channel_is_reported_without_breaking_plot_selection(qtbot):
    result = make_demo_result()
    reference = result.get_lap("baseline", "flying")
    reduced = replace(reference, channels={"speed_mps": reference.channels["speed_mps"]})
    case = replace(result.cases["baseline"], laps={"flying": reduced})
    result = replace(result, cases={"baseline": case})
    window = TelemetryViewerWindow(result)
    qtbot.addWidget(window)

    window.set_selected_channels(["control"])

    assert "control" in window.missing_channels
    assert "not available" in window.status_label.text().lower()


def test_comparison_delta_marks_slower_lap_with_positive_time(qtbot):
    window = TelemetryViewerWindow(_result_with_comparison_lap())
    qtbot.addWidget(window)

    window.set_comparison_laps("flying", "slow")

    assert window.delta_seconds == pytest.approx(0.75)
    assert "+0.750 s" in window.delta_label.text()
    assert "slower" in window.delta_label.text().lower()


def test_comparison_controls_apply_manual_offsets(qtbot):
    window = TelemetryViewerWindow(_result_with_comparison_lap())
    qtbot.addWidget(window)

    window.set_comparison_laps("flying", "slow")
    window.set_alignment_offsets(distance_offset_m=1.5, time_offset_s=0.25)

    assert window.distance_offset_spin.value() == pytest.approx(1.5)
    assert window.time_offset_spin.value() == pytest.approx(0.25)


def test_comparison_alignment_is_cached_until_inputs_change():
    calls = []

    def provider(reference, comparison):
        del reference, comparison
        calls.append(True)
        return 0.5

    model = ViewerModel(_result_with_comparison_lap(), delta_time_provider=provider)
    model.set_comparison_laps("flying", "slow")

    assert model.delta_seconds == pytest.approx(0.5)
    assert model.comparison_error is None
    assert model.comparison_status == "Comparison aligned"
    assert model.quality_summary["comparison_status"] == "Comparison aligned"
    assert len(calls) == 1

    model.set_alignment_offsets(distance_offset_m=1.0)
    assert model.delta_seconds == pytest.approx(0.5)
    assert len(calls) == 2


def test_mismatched_distance_overlap_is_reported(qtbot):
    result = make_demo_result()
    reference = result.get_lap("baseline", "flying")
    far_axis = AxisData("native", reference.axes["native"].time_s, np.array([100.0, 110.0, 120.0, 130.0]))
    far_channels = {
        channel_id: ChannelData(channel.metadata, channel.values.copy(), channel.valid.copy())
        for channel_id, channel in reference.channels.items()
    }
    far_lap = LapData("far", {"role": "comparison"}, {"native": far_axis}, far_channels, reference.track_id)
    case = replace(result.cases["baseline"], laps={"flying": reference, "far": far_lap})
    window = TelemetryViewerWindow(replace(result, cases={"baseline": case}))
    qtbot.addWidget(window)

    window.set_axis_mode("distance")
    window.set_comparison_laps("flying", "far")

    assert "comparison unavailable" in window.status_label.text().lower()
    assert "overlap" in window.status_label.text().lower()
