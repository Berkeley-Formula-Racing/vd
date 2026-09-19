from dataclasses import replace

import numpy as np
import pytest

from qss_telemetry.alignment import AlignmentError, AlignmentSettings, align_laps, calculate_delta_time
from qss_telemetry.schema import AxisData, ChannelData, ChannelMetadata, LapData


def _lap(
    lap_id,
    *,
    times=(0.0, 1.0, 2.0, 3.0),
    distances=(0.0, 10.0, 20.0, 30.0),
    values=(1.0, 1.0, 1.0, 1.0),
    valid=(True, True, True, True),
    interpolation="linear",
    track_id="track",
):
    axis = AxisData("native", np.asarray(times), np.asarray(distances))
    metadata = ChannelMetadata("signal", "Signal", "1", "native", "simulation_output", interpolation)
    channel = ChannelData(metadata, np.asarray(values, dtype=float), np.asarray(valid, dtype=bool))
    return LapData(lap_id, {}, {axis.id: axis}, {channel.metadata.id: channel}, track_id)


def test_calculate_delta_time_is_positive_when_comparison_is_slower():
    reference = _lap("reference", times=(0.0, 1.0, 2.0, 3.0))
    comparison = _lap("comparison", times=(0.0, 1.2, 2.4, 3.6))

    aligned = align_laps(reference, comparison)
    delta = calculate_delta_time(aligned)

    np.testing.assert_allclose(delta.values, [0.0, 0.2, 0.4, 0.6])
    assert delta.metadata.id == "delta_time_s"
    assert delta.metadata.sign_convention == "comparison - reference"


def test_alignment_applies_manual_distance_and_time_offsets():
    reference = _lap("reference")
    comparison = _lap(
        "comparison",
        times=(0.5, 1.5, 2.5, 3.5),
        distances=(2.0, 12.0, 22.0, 32.0),
    )

    aligned = align_laps(
        reference,
        comparison,
        AlignmentSettings(distance_offset_m=-2.0, time_offset_s=-0.5),
    )
    delta = calculate_delta_time(aligned)

    np.testing.assert_allclose(aligned.distance_m, [0.0, 10.0, 20.0, 30.0])
    np.testing.assert_allclose(delta.values, [0.0, 0.0, 0.0, 0.0])


def test_alignment_preserves_invalid_runs_and_uses_previous_for_discrete_channels():
    reference = _lap(
        "reference",
        values=(1.0, np.nan, np.nan, 4.0),
        valid=(True, False, False, True),
        interpolation="linear",
    )
    comparison = _lap(
        "comparison",
        values=(1.0, 2.0, 3.0, 4.0),
        interpolation="previous",
    )

    aligned = align_laps(reference, comparison)

    assert aligned.reference_channels["signal"].valid.tolist() == [True, False, False, True]
    assert aligned.comparison_channels["signal"].valid.tolist() == [True, True, True, True]
    np.testing.assert_allclose(aligned.comparison_channels["signal"].values, [1.0, 2.0, 3.0, 4.0])


def test_alignment_rejects_mismatched_tracks():
    reference = _lap("reference", track_id="track-a")
    comparison = _lap("comparison", track_id="track-b")

    with pytest.raises(AlignmentError, match="track"):
        align_laps(reference, comparison)
