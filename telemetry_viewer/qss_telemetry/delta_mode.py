"""Helpers for rendering telemetry channels relative to a datum lap."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .schema import LapData


@dataclass(frozen=True, slots=True)
class DeltaSeries:
    """A common x-axis and active-minus-datum waveform."""

    x: np.ndarray
    values: np.ndarray
    valid: np.ndarray


def _axis_values(lap: LapData, channel_id: str, axis_mode: str) -> np.ndarray:
    channel = lap.channels[channel_id]
    axis = lap.axes[channel.metadata.axis_id]
    return np.asarray(axis.distance_m if axis_mode == "distance" else axis.time_s, dtype=float)


def compute_delta_series(
    active_lap: LapData,
    datum_lap: LapData,
    channel_id: str,
    *,
    axis_mode: str = "time",
) -> DeltaSeries:
    """Return ``active - datum`` on the union of both channel axes.

    Both channels are linearly interpolated onto the union of their sampled
    coordinates.  Samples outside the valid overlap of the two channels are
    represented as NaN, while the viewer's cursor remains tied to the active
    lap.  Using the union keeps a slower or differently sampled datum from
    disappearing between active-lap samples.
    """

    if channel_id not in active_lap.channels:
        raise KeyError(f"active lap is missing channel: {channel_id}")
    if channel_id not in datum_lap.channels:
        raise KeyError(f"datum lap is missing channel: {channel_id}")

    active_channel = active_lap.channels[channel_id]
    datum_channel = datum_lap.channels[channel_id]
    if active_channel.metadata.unit != datum_channel.metadata.unit:
        raise ValueError(
            f"channel {channel_id} has incompatible units: "
            f"{active_channel.metadata.unit!r} vs {datum_channel.metadata.unit!r}"
        )

    active_x = _axis_values(active_lap, channel_id, axis_mode)
    datum_x = _axis_values(datum_lap, channel_id, axis_mode)
    active_values = np.asarray(active_channel.values, dtype=float)
    datum_values = np.asarray(datum_channel.values, dtype=float)
    active_valid = np.asarray(active_channel.valid, dtype=bool) & np.isfinite(active_values)
    datum_valid = np.asarray(datum_channel.valid, dtype=bool) & np.isfinite(datum_values) & np.isfinite(datum_x)

    finite_active_x = active_x[np.isfinite(active_x)]
    finite_datum_x = datum_x[np.isfinite(datum_x)]
    common_x = np.unique(np.concatenate((finite_active_x, finite_datum_x)))
    delta = np.full(common_x.shape, np.nan, dtype=float)
    valid = np.zeros(common_x.shape, dtype=bool)
    if not np.any(active_valid) or not np.any(datum_valid) or common_x.size == 0:
        return DeltaSeries(common_x, delta, valid)

    def _interpolate_valid(
        x_values: np.ndarray,
        values: np.ndarray,
        mask: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        interpolated = np.full(common_x.shape, np.nan, dtype=float)
        supported = np.zeros(common_x.shape, dtype=bool)
        valid_indices = np.flatnonzero(mask)
        runs = np.split(valid_indices, np.flatnonzero(np.diff(valid_indices) > 1) + 1)
        for run in runs:
            if run.size == 0:
                continue
            run_x, unique_indices = np.unique(x_values[run], return_index=True)
            run_values = values[run][unique_indices]
            if run_x.size == 1:
                exact = common_x == run_x[0]
                interpolated[exact] = run_values[0]
                supported[exact] = True
                continue
            inside = (common_x >= run_x[0]) & (common_x <= run_x[-1])
            interpolated[inside] = np.interp(common_x[inside], run_x, run_values)
            supported[inside] = True
        return interpolated, supported

    active_interpolated, active_supported = _interpolate_valid(active_x, active_values, active_valid)
    datum_interpolated, datum_supported = _interpolate_valid(datum_x, datum_values, datum_valid)
    valid = active_supported & datum_supported
    if np.any(valid):
        delta[valid] = active_interpolated[valid] - datum_interpolated[valid]
    return DeltaSeries(common_x, delta, valid)
