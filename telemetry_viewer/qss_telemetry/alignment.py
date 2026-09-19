"""Distance/time alignment for v1 QSS telemetry laps.

Alignment is deliberately local and deterministic.  It uses a common
distance grid by default, restricts that grid to the physical overlap, and
resamples each channel according to its metadata.  Invalid source runs stay
invalid in the result, so interpolation never invents samples across a
telemetry gap.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Literal, Mapping, Sequence

import numpy as np

from .schema import AxisData, ChannelData, ChannelMetadata, LapData, SchemaError


class AlignmentError(ValueError):
    """Raised when two laps cannot be aligned under v1 rules."""


@dataclass(frozen=True, slots=True)
class AlignmentSettings:
    """Options controlling common-coordinate alignment.

    ``mode`` is ``"distance"`` by default and may be ``"time"`` for a
    time-based comparison.  Generic ``distance_offset_m`` and
    ``time_offset_s`` apply to the comparison side; explicit side-specific
    fields are available when both laps need independent corrections.  A
    positive comparison offset shifts that lap's coordinates later/right in
    the common frame.  The optional grid/step fields are conveniences for a
    plotter; absent all of them, the sorted union of source coordinates in
    the overlap is used.
    """

    mode: Literal["distance", "time"] = "distance"
    distance_offset_m: float = 0.0
    time_offset_s: float = 0.0
    reference_distance_offset_m: float = 0.0
    comparison_distance_offset_m: float | None = None
    reference_time_offset_s: float = 0.0
    comparison_time_offset_s: float | None = None
    grid: Sequence[float] | None = None
    distance_grid: Sequence[float] | None = None
    time_grid: Sequence[float] | None = None
    step: float | None = None
    distance_step_m: float | None = None
    time_step_s: float | None = None
    sample_count: int | None = None
    axis_id: str | None = None
    reference_axis_id: str | None = None
    comparison_axis_id: str | None = None
    # Readable aliases retained for callers that use these names.
    coordinate: str | None = None
    coordinate_mode: str | None = None
    basis: str | None = None
    manual_distance_offset_m: float = 0.0
    manual_time_offset_s: float = 0.0

    def __post_init__(self) -> None:
        mode = self.effective_mode
        if mode not in {"distance", "time"}:
            raise ValueError("alignment mode must be 'distance' or 'time'")
        for name in (
            "distance_offset_m",
            "time_offset_s",
            "reference_distance_offset_m",
            "reference_time_offset_s",
            "manual_distance_offset_m",
            "manual_time_offset_s",
        ):
            if not np.isfinite(float(getattr(self, name))):
                raise ValueError(f"{name} must be finite")
        for name in ("comparison_distance_offset_m", "comparison_time_offset_s"):
            value = getattr(self, name)
            if value is not None and not np.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
        if self.sample_count is not None and (isinstance(self.sample_count, bool) or self.sample_count < 1):
            raise ValueError("sample_count must be a positive integer")
        for name in ("step", "distance_step_m", "time_step_s"):
            value = getattr(self, name)
            if value is not None and (not np.isfinite(float(value)) or float(value) <= 0):
                raise ValueError(f"{name} must be a positive finite number")
        for name in ("grid", "distance_grid", "time_grid"):
            value = getattr(self, name)
            if value is not None:
                array = np.asarray(value, dtype=float)
                if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
                    raise ValueError(f"{name} must be a non-empty finite one-dimensional sequence")

    @property
    def effective_mode(self) -> str:
        raw = self.coordinate_mode or self.coordinate or self.basis or self.mode
        raw = str(raw).lower()
        if raw in {"distance", "distance_m", "physical_distance", "s"}:
            return "distance"
        if raw in {"time", "time_s", "seconds", "t"}:
            return "time"
        return raw

    @property
    def reference_distance_offset(self) -> float:
        return float(self.reference_distance_offset_m) - 0.0

    @property
    def comparison_distance_offset(self) -> float:
        generic = self.distance_offset_m if self.comparison_distance_offset_m is None else self.comparison_distance_offset_m
        return float(generic) + float(self.manual_distance_offset_m)

    @property
    def reference_time_offset(self) -> float:
        return float(self.reference_time_offset_s)

    @property
    def comparison_time_offset(self) -> float:
        generic = self.time_offset_s if self.comparison_time_offset_s is None else self.comparison_time_offset_s
        return float(generic) + float(self.manual_time_offset_s)

    def grid_for_mode(self) -> Sequence[float] | None:
        if self.effective_mode == "distance":
            return self.distance_grid if self.distance_grid is not None else self.grid
        return self.time_grid if self.time_grid is not None else self.grid

    def step_for_mode(self) -> float | None:
        if self.effective_mode == "distance":
            return self.distance_step_m if self.distance_step_m is not None else self.step
        return self.time_step_s if self.time_step_s is not None else self.step


@dataclass(frozen=True, slots=True)
class AlignedPair:
    """Two laps sampled on one common coordinate grid.

    ``distance_m`` is the canonical physical coordinate, while
    ``reference_time_s`` and ``comparison_time_s`` retain each lap's sampled
    time.  Channel mappings contain resampled copies keyed by channel id;
    channels absent from one source are therefore simply absent from that
    side's mapping.  ``axis`` is a synthetic ``alignment`` axis suitable for
    plotting those copies.
    """

    reference: LapData
    comparison: LapData
    axis: AxisData
    distance_m: np.ndarray
    reference_time_s: np.ndarray
    comparison_time_s: np.ndarray
    reference_channels: Mapping[str, ChannelData]
    comparison_channels: Mapping[str, ChannelData]
    settings: AlignmentSettings

    @property
    def common_axis(self) -> AxisData:
        return self.axis

    @property
    def time_s(self) -> np.ndarray:
        return self.axis.time_s

    @property
    def coordinate(self) -> np.ndarray:
        return self.distance_m if self.settings.effective_mode == "distance" else self.axis.time_s

    @property
    def coordinate_mode(self) -> str:
        return self.settings.effective_mode

    @property
    def ref(self) -> LapData:
        return self.reference

    @property
    def comp(self) -> LapData:
        return self.comparison

    @property
    def channels(self) -> Mapping[str, tuple[ChannelData | None, ChannelData | None]]:
        ids = set(self.reference_channels) | set(self.comparison_channels)
        return {channel_id: (self.reference_channels.get(channel_id), self.comparison_channels.get(channel_id)) for channel_id in ids}


def _coerce_settings(settings: AlignmentSettings | Mapping[str, Any] | None) -> AlignmentSettings:
    if settings is None:
        return AlignmentSettings()
    if isinstance(settings, AlignmentSettings):
        return settings
    if isinstance(settings, Mapping):
        try:
            return AlignmentSettings(**dict(settings))
        except TypeError as error:
            raise AlignmentError(f"invalid alignment settings: {error}") from error
    raise TypeError("settings must be AlignmentSettings, a mapping, or None")


def _pick_axis(lap: LapData, settings: AlignmentSettings, side: Literal["reference", "comparison"]) -> AxisData:
    requested = settings.reference_axis_id if side == "reference" else settings.comparison_axis_id
    requested = requested or settings.axis_id
    if requested is not None:
        try:
            return lap.axes[requested]
        except KeyError as error:
            raise AlignmentError(f"{side} lap has no axis {requested!r}") from error
    if "native" in lap.axes:
        return lap.axes["native"]
    return next(iter(lap.axes.values()))


def _compressed_pairs(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Collapse duplicate monotonic coordinates, keeping the latest value."""

    if x.size == 0:
        return x, y
    keep_x: list[float] = [float(x[0])]
    keep_y: list[float] = [float(y[0])]
    for coordinate, value in zip(x[1:], y[1:]):
        coordinate = float(coordinate)
        if coordinate == keep_x[-1]:
            keep_y[-1] = float(value)
        else:
            keep_x.append(coordinate)
            keep_y.append(float(value))
    return np.asarray(keep_x, dtype=float), np.asarray(keep_y, dtype=float)


def _axis_interpolate(source_x: np.ndarray, source_y: np.ndarray, target_x: np.ndarray) -> np.ndarray:
    x, y = _compressed_pairs(np.asarray(source_x, dtype=float), np.asarray(source_y, dtype=float))
    if x.size == 0:
        return np.full(target_x.shape, np.nan, dtype=float)
    if x.size == 1:
        return np.full(target_x.shape, y[0], dtype=float)
    return np.interp(target_x, x, y)


def _valid_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(indices) > 1)
    starts = np.r_[0, breaks + 1]
    ends = np.r_[breaks, indices.size - 1]
    return [(int(indices[start]), int(indices[end])) for start, end in zip(starts, ends)]


def _resample_channel(channel: ChannelData, source_coordinate: np.ndarray, target: np.ndarray) -> ChannelData:
    source_coordinate = np.asarray(source_coordinate, dtype=float)
    values = np.asarray(channel.values, dtype=float)
    source_valid = np.asarray(channel.valid, dtype=bool)
    if source_coordinate.ndim != 1 or source_coordinate.size != values.size:
        raise AlignmentError(f"channel {channel.metadata.id!r} coordinate length does not match values")
    finite_source = np.isfinite(source_coordinate) & np.isfinite(values)
    mask = source_valid & finite_source
    output = np.full(target.shape, np.nan, dtype=float)
    output_valid = np.zeros(target.shape, dtype=bool)
    mode = channel.metadata.interpolation
    for start, end in _valid_runs(mask):
        x = source_coordinate[start : end + 1]
        y = values[start : end + 1]
        # A monotonic axis can repeat a coordinate.  Keeping the latest value
        # gives deterministic exact/previous behavior and avoids zero-width
        # linear interpolation intervals.
        x, y = _compressed_pairs(x, y)
        if x.size == 0:
            continue
        if mode == "linear":
            if x.size > 1:
                selected = (target >= x[0]) & (target <= x[-1])
                if np.any(selected):
                    output[selected] = np.interp(target[selected], x, y)
                    output_valid[selected] = True
            else:
                selected = np.isclose(target, x[0], rtol=1e-10, atol=1e-12)
                output[selected] = y[0]
                output_valid[selected] = True
        elif mode == "previous" and x.size >= 1:
            selected = (target >= x[0]) & (target <= x[-1])
            if np.any(selected):
                positions = np.searchsorted(x, target[selected], side="right") - 1
                usable = positions >= 0
                selected_indices = np.flatnonzero(selected)
                if np.any(usable):
                    output[selected_indices[usable]] = y[positions[usable]]
                    output_valid[selected_indices[usable]] = True
        elif mode == "none":
            for coordinate, value in zip(x, y):
                selected = np.isclose(target, coordinate, rtol=1e-10, atol=1e-12)
                output[selected] = value
                output_valid[selected] = True
        else:
            raise AlignmentError(f"unsupported interpolation mode {mode!r} for channel {channel.metadata.id!r}")
    metadata = replace(channel.metadata, axis_id="alignment")
    try:
        return ChannelData(metadata, output, output_valid)
    except SchemaError as error:
        raise AlignmentError(f"could not resample channel {channel.metadata.id!r}: {error}") from error


def _common_grid(reference_coordinate: np.ndarray, comparison_coordinate: np.ndarray, settings: AlignmentSettings) -> np.ndarray:
    ref = np.asarray(reference_coordinate, dtype=float)
    comp = np.asarray(comparison_coordinate, dtype=float)
    if ref.size == 0 or comp.size == 0:
        raise AlignmentError("cannot align an empty axis")
    start = max(float(ref[0]), float(comp[0]))
    end = min(float(ref[-1]), float(comp[-1]))
    if end < start:
        raise AlignmentError("laps have no overlapping alignment interval")
    grid = settings.grid_for_mode()
    if grid is not None:
        target = np.asarray(grid, dtype=float)
        target = target[(target >= start) & (target <= end)]
        if target.size == 0:
            raise AlignmentError("requested alignment grid has no points in the overlapping interval")
    else:
        step = settings.step_for_mode()
        if step is not None:
            count = int(np.floor((end - start) / step)) + 1
            target = start + np.arange(max(count, 1), dtype=float) * step
            if target[-1] < end:
                target = np.r_[target, end]
            else:
                target[-1] = end
        elif settings.sample_count is not None:
            target = np.linspace(start, end, int(settings.sample_count), dtype=float)
        else:
            target = np.r_[ref[(ref >= start) & (ref <= end)], comp[(comp >= start) & (comp <= end)], start, end]
    target = np.unique(target)
    if target.size == 0:
        raise AlignmentError("laps have no usable alignment samples")
    return target


def _lap_coordinate(axis: AxisData, settings: AlignmentSettings, side: Literal["reference", "comparison"]) -> np.ndarray:
    if settings.effective_mode == "distance":
        offset = settings.reference_distance_offset if side == "reference" else settings.comparison_distance_offset
        return np.asarray(axis.distance_m, dtype=float) + offset
    offset = settings.reference_time_offset if side == "reference" else settings.comparison_time_offset
    return np.asarray(axis.time_s, dtype=float) + offset


def _channel_coordinate(axis: AxisData, settings: AlignmentSettings, side: Literal["reference", "comparison"]) -> np.ndarray:
    return _lap_coordinate(axis, settings, side)


def align_laps(
    reference: LapData,
    comparison: LapData,
    settings: AlignmentSettings | Mapping[str, Any] | None = None,
) -> AlignedPair:
    """Resample two laps onto the overlapping v1 distance or time interval."""

    if not isinstance(reference, LapData) or not isinstance(comparison, LapData):
        raise TypeError("reference and comparison must be LapData")
    options = _coerce_settings(settings)
    if reference.track_id != comparison.track_id:
        raise AlignmentError(f"cannot align laps from different tracks: {reference.track_id!r} and {comparison.track_id!r}")
    reference_axis = _pick_axis(reference, options, "reference")
    comparison_axis = _pick_axis(comparison, options, "comparison")
    reference_coordinate = _lap_coordinate(reference_axis, options, "reference")
    comparison_coordinate = _lap_coordinate(comparison_axis, options, "comparison")
    target = _common_grid(reference_coordinate, comparison_coordinate, options)

    if options.effective_mode == "distance":
        distance_m = target
        reference_time = _axis_interpolate(reference_coordinate, reference_axis.time_s + options.reference_time_offset, target)
        comparison_time = _axis_interpolate(comparison_coordinate, comparison_axis.time_s + options.comparison_time_offset, target)
        common_time = reference_time
    else:
        reference_time = _axis_interpolate(reference_coordinate, reference_axis.time_s + options.reference_time_offset, target)
        comparison_time = _axis_interpolate(comparison_coordinate, comparison_axis.time_s + options.comparison_time_offset, target)
        distance_source = np.asarray(reference_axis.distance_m, dtype=float) + options.reference_distance_offset
        distance_m = _axis_interpolate(reference_coordinate, distance_source, target)
        common_time = target

    aligned_axis = AxisData("alignment", common_time, distance_m)
    reference_channels: dict[str, ChannelData] = {}
    comparison_channels: dict[str, ChannelData] = {}
    for channel_id, channel in reference.channels.items():
        source_axis = reference.axes[channel.metadata.axis_id]
        reference_channels[channel_id] = _resample_channel(channel, _channel_coordinate(source_axis, options, "reference"), target)
    for channel_id, channel in comparison.channels.items():
        source_axis = comparison.axes[channel.metadata.axis_id]
        comparison_channels[channel_id] = _resample_channel(channel, _channel_coordinate(source_axis, options, "comparison"), target)
    return AlignedPair(
        reference=reference,
        comparison=comparison,
        axis=aligned_axis,
        distance_m=distance_m,
        reference_time_s=reference_time,
        comparison_time_s=comparison_time,
        reference_channels=reference_channels,
        comparison_channels=comparison_channels,
        settings=options,
    )


def calculate_delta_time(aligned_pair: AlignedPair) -> ChannelData:
    """Return ``comparison_time - reference_time`` on an aligned grid.

    Positive values mean the comparison lap is slower at that coordinate.
    ``NaN`` and ``valid=False`` are retained wherever either aligned time is
    unavailable.
    """

    if not isinstance(aligned_pair, AlignedPair):
        raise TypeError("aligned_pair must be an AlignedPair")
    reference_time = np.asarray(aligned_pair.reference_time_s, dtype=float)
    comparison_time = np.asarray(aligned_pair.comparison_time_s, dtype=float)
    if reference_time.shape != comparison_time.shape or reference_time.shape != aligned_pair.distance_m.shape:
        raise AlignmentError("aligned time and distance arrays must have equal shapes")
    valid = np.isfinite(reference_time) & np.isfinite(comparison_time)
    values = np.full(reference_time.shape, np.nan, dtype=float)
    values[valid] = comparison_time[valid] - reference_time[valid]
    metadata = ChannelMetadata(
        id="delta_time_s",
        label="Delta time",
        unit="s",
        axis_id=aligned_pair.axis.id,
        origin="derived",
        interpolation="linear",
        description="Comparison time minus reference time at the aligned coordinate.",
        coordinate_frame="vehicle",
        sign_convention="comparison - reference",
    )
    return ChannelData(metadata, values, valid)


__all__ = [
    "AlignedPair",
    "AlignmentError",
    "AlignmentSettings",
    "align_laps",
    "calculate_delta_time",
]
