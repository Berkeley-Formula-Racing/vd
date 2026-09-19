from __future__ import annotations

import numpy as np
import pytest

from qss_telemetry.delta_mode import compute_delta_series
from qss_telemetry.fixtures import make_demo_result
from qss_telemetry.schema import AxisData, ChannelData, LapData


def _lap_with_speed(axis: AxisData, values: list[float], valid: list[bool] | None = None) -> LapData:
    base = make_demo_result().get_lap("baseline", "flying")
    speed = base.channels["speed_mps"]
    mask = np.ones(len(values), dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
    channel = ChannelData(speed.metadata, np.where(mask, values, np.nan), mask)
    return LapData("lap", {}, {"native": axis}, {"speed_mps": channel}, base.track_id)


def test_delta_is_interpolated_on_the_active_axis_and_invalid_overlap_is_nan():
    active_axis = AxisData("native", np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.0, 10.0, 20.0, 30.0]))
    datum_axis = AxisData("native", np.array([0.5, 1.5, 2.5]), np.array([5.0, 15.0, 25.0]))
    active = _lap_with_speed(active_axis, [12.0, 14.0, 16.0, 18.0])
    datum = _lap_with_speed(datum_axis, [10.0, 12.0, 14.0])

    result = compute_delta_series(active, datum, "speed_mps")

    assert result.x.tolist() == [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    assert result.values.tolist() == pytest.approx([np.nan, 3.0, 3.0, 3.0, 3.0, 3.0, np.nan], nan_ok=True)
    assert result.valid.tolist() == [False, True, True, True, True, True, False]


def test_delta_rejects_unit_mismatch():
    result = make_demo_result()
    active = result.get_lap("baseline", "flying")
    speed = active.channels["speed_mps"]
    datum_channel = ChannelData(speed.metadata.__class__(**{**speed.metadata.as_dict(), "unit": "mph"}), speed.values, speed.valid)
    datum = LapData("datum", {}, active.axes, {"speed_mps": datum_channel}, active.track_id)

    with pytest.raises(ValueError, match="incompatible units"):
        compute_delta_series(active, datum, "speed_mps")
