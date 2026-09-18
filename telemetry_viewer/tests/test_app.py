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

from qss_telemetry.app import TelemetryViewerWindow
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


def test_fixture_result_is_visible_in_browser_and_default_channels(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)

    assert window.case_lap_tree.topLevelItemCount() == 1
    assert window.case_lap_tree.topLevelItem(0).text(0) == "Baseline"
    assert window.channel_selector.is_checked("speed_mps")
    assert window.channel_selector.is_checked("long_accel_mps2")
    assert window.channel_selector.is_checked("lat_accel_mps2")
    assert window.channel_selector.is_checked("gear")


def test_cursor_updates_all_plots_and_schematic_map_marker(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)

    window.set_cursor_index(2)

    assert window.cursor_index == 2
    assert window.map_widget.marker_distance_m == pytest.approx(20.0)
    assert window.map_widget.quality_label.text() == "Schematic map"
    assert all(line.value() == pytest.approx(2.0) for line in window.plot_stack.cursor_lines)


def test_playback_step_advances_cursor_and_map_deterministically(qtbot):
    window = TelemetryViewerWindow(make_demo_result())
    qtbot.addWidget(window)
    window.set_cursor_index(0)

    window.playback_controller.step()

    assert window.cursor_index == 1
    assert window.map_widget.marker_distance_m == pytest.approx(10.0)


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
