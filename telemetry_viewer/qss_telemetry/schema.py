"""Stable, dependency-light data contracts for QSS telemetry version 1.0."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping

import numpy as np
from numpy.typing import NDArray

FORMAT_NAME = "qss-telemetry"
SCHEMA_VERSION = (1, 0)
MANIFEST_DATASET = "/manifest_json"

Origin = Literal["simulation_output", "derived", "qss_reconstructed", "measured"]
Interpolation = Literal["linear", "previous", "none"]


class SchemaError(ValueError):
    """Raised when an in-memory result violates the v1 telemetry contract."""


def _as_float_vector(name: str, values: NDArray[np.floating[Any]] | np.ndarray) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise SchemaError(f"{name} must be one-dimensional")
    return array


@dataclass(frozen=True, slots=True)
class ChannelMetadata:
    """Self-describing metadata persisted as a UTF-8 JSON object."""

    id: str
    label: str
    unit: str
    axis_id: str
    origin: Origin
    interpolation: Interpolation
    description: str = ""
    coordinate_frame: str = "vehicle"
    sign_convention: str = ""

    def as_dict(self) -> dict[str, str]:
        return {
            "id": self.id,
            "label": self.label,
            "unit": self.unit,
            "axis_id": self.axis_id,
            "origin": self.origin,
            "interpolation": self.interpolation,
            "description": self.description,
            "coordinate_frame": self.coordinate_frame,
            "sign_convention": self.sign_convention,
        }


@dataclass(frozen=True, slots=True)
class AxisData:
    """A monotonic time-distance coordinate pair for a sampled lap."""

    id: str
    time_s: NDArray[np.float64]
    distance_m: NDArray[np.float64]

    def __post_init__(self) -> None:
        time_s = _as_float_vector("time_s", self.time_s)
        distance_m = _as_float_vector("distance_m", self.distance_m)
        if time_s.size == 0 or time_s.size != distance_m.size:
            raise SchemaError("an axis needs equal-length, non-empty time_s and distance_m")
        if not (np.all(np.isfinite(time_s)) and np.all(np.isfinite(distance_m))):
            raise SchemaError("axis coordinates must be finite")
        if np.any(np.diff(time_s) < 0) or np.any(np.diff(distance_m) < 0):
            raise SchemaError("axis coordinates must be monotonic")
        object.__setattr__(self, "time_s", time_s)
        object.__setattr__(self, "distance_m", distance_m)


@dataclass(frozen=True, slots=True)
class ChannelData:
    """One sampled channel, including explicit validity rather than zero filling."""

    metadata: ChannelMetadata
    values: NDArray[np.float64]
    valid: NDArray[np.bool_]

    def __post_init__(self) -> None:
        values = _as_float_vector("channel values", self.values)
        valid = np.asarray(self.valid, dtype=bool)
        if valid.ndim != 1 or values.shape != valid.shape:
            raise SchemaError("channel values and valid mask must be equal-length vectors")
        if np.any(~valid & ~np.isnan(values)):
            raise SchemaError("invalid samples must be NaN")
        if np.any(valid & ~np.isfinite(values)):
            raise SchemaError("valid samples must be finite")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "valid", valid)


@dataclass(frozen=True, slots=True)
class TrackData:
    id: str
    metadata: Mapping[str, Any]
    distance_m: NDArray[np.float64]
    curvature_per_m: NDArray[np.float64]
    x_m: NDArray[np.float64] | None = None
    y_m: NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        distance_m = _as_float_vector("track distance_m", self.distance_m)
        curvature = _as_float_vector("track curvature_per_m", self.curvature_per_m)
        if distance_m.size == 0 or distance_m.size != curvature.size:
            raise SchemaError("track distance and curvature must be equal-length and non-empty")
        if np.any(np.diff(distance_m) < 0):
            raise SchemaError("track distance must be monotonic")
        for name in ("x_m", "y_m"):
            candidate = getattr(self, name)
            if candidate is not None:
                candidate = _as_float_vector(f"track {name}", candidate)
                if candidate.shape != distance_m.shape:
                    raise SchemaError(f"track {name} must match track distance")
                object.__setattr__(self, name, candidate)
        if (self.x_m is None) != (self.y_m is None):
            raise SchemaError("track x_m and y_m must either both exist or both be absent")
        object.__setattr__(self, "distance_m", distance_m)
        object.__setattr__(self, "curvature_per_m", curvature)


@dataclass(frozen=True, slots=True)
class LapData:
    id: str
    metadata: Mapping[str, Any]
    axes: Mapping[str, AxisData]
    channels: Mapping[str, ChannelData]
    track_id: str

    def __post_init__(self) -> None:
        if not self.axes:
            raise SchemaError("a lap must contain at least one axis")
        for channel_id, channel in self.channels.items():
            if channel_id != channel.metadata.id:
                raise SchemaError("channel map key must match metadata id")
            axis = self.axes.get(channel.metadata.axis_id)
            if axis is None:
                raise SchemaError(f"channel {channel_id} refers to a missing axis")
            if channel.values.size != axis.time_s.size:
                raise SchemaError(f"channel {channel_id} length does not match its axis")


@dataclass(frozen=True, slots=True)
class CaseData:
    id: str
    metadata: Mapping[str, Any]
    setup: Mapping[str, Any]
    laps: Mapping[str, LapData]
    envelope: Mapping[str, ChannelData] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class ResultFile:
    manifest: Mapping[str, Any]
    tracks: Mapping[str, TrackData]
    cases: Mapping[str, CaseData]
    source_path: str | None = None

    def get_lap(self, case_id: str, lap_id: str) -> LapData:
        try:
            return self.cases[case_id].laps[lap_id]
        except KeyError as error:
            raise KeyError(f"unknown case/lap: {case_id}/{lap_id}") from error
