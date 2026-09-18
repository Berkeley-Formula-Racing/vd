"""Validation helpers for QSS telemetry files and in-memory laps.

The immutable classes in :mod:`qss_telemetry.schema` perform the final
in-memory checks.  This module adds non-throwing diagnostics that can be
shown by the viewer before a file is loaded.  A report is valid when it has
no error-severity issues; warnings are retained for useful, non-fatal
omissions such as optional manifest provenance fields.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping

import h5py
import numpy as np

from .schema import FORMAT_NAME, SCHEMA_VERSION, AxisData, CaseData, ChannelData, LapData, ResultFile, TrackData


Severity = Literal["error", "warning"]


@dataclass(frozen=True, slots=True)
class ValidationIssue:
    """One compact, UI-friendly validation diagnostic."""

    path: str
    message: str
    code: str = "invalid"
    severity: Severity = "error"

    def __str__(self) -> str:
        location = f"{self.path}: " if self.path else ""
        return f"{location}{self.message}"

    def __contains__(self, value: object) -> bool:
        """Allow convenient ``"text" in issue`` assertions in UI/tests."""

        return str(value) in str(self)


@dataclass(frozen=True, slots=True)
class ValidationReport:
    """Collection of validation diagnostics.

    ``issues`` is intentionally a tuple so reports can safely cross the
    worker/UI boundary.  ``valid`` and its aliases make the report pleasant
    to consume from either Python or a Qt model.
    """

    issues: tuple[ValidationIssue, ...] = ()

    @property
    def valid(self) -> bool:
        return not any(issue.severity == "error" for issue in self.issues)

    @property
    def ok(self) -> bool:
        return self.valid

    @property
    def is_valid(self) -> bool:
        return self.valid

    def __bool__(self) -> bool:
        return self.valid

    @property
    def errors(self) -> tuple[ValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.severity == "error")

    @property
    def warnings(self) -> tuple[ValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.severity == "warning")

    @property
    def error_messages(self) -> tuple[str, ...]:
        return tuple(str(issue) for issue in self.errors)

    @property
    def warning_messages(self) -> tuple[str, ...]:
        return tuple(str(issue) for issue in self.warnings)

    @property
    def summary(self) -> str:
        if self.valid:
            return "valid" if not self.issues else f"valid with {len(self.warnings)} warning(s)"
        return "; ".join(str(issue) for issue in self.errors)


def _report(issues: Iterable[ValidationIssue]) -> ValidationReport:
    return ValidationReport(tuple(issues))


def _issue(path: str, message: str, code: str = "invalid", severity: Severity = "error") -> ValidationIssue:
    return ValidationIssue(path=path, message=message, code=code, severity=severity)


def _check_json_dataset(dataset: h5py.Dataset, path: str, issues: list[ValidationIssue], *, object_name: str) -> Any | None:
    """Decode a JSON dataset for validation without letting malformed bytes escape."""

    array: np.ndarray | None = None
    raw: Any
    try:
        raw = dataset[()]
        array = np.asarray(raw)
    except Exception as error:  # h5py can raise several backend-specific errors
        issues.append(_issue(path, f"could not read {object_name} JSON: {error}", "json_read"))
        return None

    if array.ndim != 1 or array.dtype != np.dtype(np.uint8):
        issues.append(_issue(path, f"{object_name} JSON must be a one-dimensional uint8 dataset", "json_encoding"))
        return None
    try:
        value = bytes(array.tolist()).decode("utf-8")
        value = __import__("json").loads(value)
    except (UnicodeDecodeError, ValueError, TypeError) as error:
        issues.append(_issue(path, f"malformed {object_name} JSON: {error}", "json_decode"))
        return None
    if not isinstance(value, Mapping):
        issues.append(_issue(path, f"{object_name} JSON must contain an object", "json_type"))
    return value


def _check_dataset(group: h5py.Group, name: str, path: str, issues: list[ValidationIssue]) -> h5py.Dataset | None:
    if name not in group:
        issues.append(_issue(f"{path}/{name}", "required dataset is missing", "missing_dataset"))
        return None
    item = group[name]
    if not isinstance(item, h5py.Dataset):
        issues.append(_issue(f"{path}/{name}", "expected a dataset", "dataset_type"))
        return None
    return item


def _read_array(dataset: h5py.Dataset, path: str, issues: list[ValidationIssue]) -> np.ndarray | None:
    try:
        return np.asarray(dataset[()])
    except Exception as error:
        issues.append(_issue(path, f"could not read dataset: {error}", "dataset_read"))
        return None


def _check_numeric_vector(dataset: h5py.Dataset, path: str, issues: list[ValidationIssue], *, name: str) -> np.ndarray | None:
    array = _read_array(dataset, path, issues)
    if array is None:
        return None
    if array.ndim != 1:
        issues.append(_issue(path, f"{name} must be one-dimensional", "shape"))
        return array
    if not np.issubdtype(array.dtype, np.number):
        issues.append(_issue(path, f"{name} must be numeric", "dtype"))
        return array
    try:
        converted = np.asarray(array, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        issues.append(_issue(path, f"{name} must be float-compatible: {error}", "dtype"))
        return array
    if not np.all(np.isfinite(converted)):
        issues.append(_issue(path, f"{name} must contain only finite values", "finite"))
    if converted.size == 0:
        issues.append(_issue(path, f"{name} must be non-empty", "empty"))
    if converted.size > 1 and np.any(np.diff(converted) < 0):
        issues.append(_issue(path, f"{name} must be monotonic non-decreasing", "monotonic"))
    return converted


def _check_values_and_mask(
    values_dataset: h5py.Dataset | None,
    valid_dataset: h5py.Dataset | None,
    values_path: str,
    valid_path: str,
    issues: list[ValidationIssue],
) -> tuple[np.ndarray | None, np.ndarray | None]:
    values = _read_array(values_dataset, values_path, issues) if values_dataset is not None else None
    valid = _read_array(valid_dataset, valid_path, issues) if valid_dataset is not None else None
    if values is not None:
        if values.ndim != 1:
            issues.append(_issue(values_path, "channel values must be one-dimensional", "shape"))
        elif not np.issubdtype(values.dtype, np.number):
            issues.append(_issue(values_path, "channel values must be numeric", "dtype"))
        else:
            values = np.asarray(values, dtype=np.float64)
    if valid is not None:
        if valid.ndim != 1:
            issues.append(_issue(valid_path, "valid mask must be one-dimensional", "shape"))
        if valid.dtype != np.dtype(np.uint8):
            issues.append(_issue(valid_path, "valid mask must use uint8 values", "dtype"))
        elif np.any((valid != 0) & (valid != 1)):
            issues.append(_issue(valid_path, "valid mask values must be 0 or 1", "mask_value"))
        valid = np.asarray(valid, dtype=bool)
    if values is not None and valid is not None and values.ndim == valid.ndim == 1:
        if values.shape != valid.shape:
            issues.append(_issue(values_path, "channel values and valid mask must have equal shapes", "shape"))
        else:
            invalid = ~valid
            if np.any(invalid & ~np.isnan(values)):
                issues.append(_issue(values_path, "invalid samples must be NaN", "nan_mask"))
            if np.any(valid & ~np.isfinite(values)):
                issues.append(_issue(values_path, "valid samples must be finite", "finite"))
    return values, valid


def _metadata_fields(metadata: Mapping[str, Any] | None, path: str, issues: list[ValidationIssue], *, channel_id: str | None = None) -> None:
    required = ("id", "label", "unit", "axis_id", "origin", "interpolation", "description", "coordinate_frame", "sign_convention")
    if metadata is None:
        return
    for key in required:
        if key not in metadata:
            issues.append(_issue(f"{path}/{key}", "required metadata field is missing", "missing_metadata"))
        elif not isinstance(metadata[key], str):
            issues.append(_issue(f"{path}/{key}", "channel metadata fields must be strings", "metadata_type"))
        elif key in {"id", "axis_id", "origin", "interpolation"} and not metadata[key]:
            issues.append(_issue(f"{path}/{key}", "channel metadata field must not be empty", "metadata_value"))
    if channel_id is not None and metadata.get("id") != channel_id:
        issues.append(_issue(path, f"channel metadata id must match group id ({channel_id!r})", "metadata_id"))
    if "origin" in metadata and metadata["origin"] not in {"simulation_output", "derived", "qss_reconstructed", "measured"}:
        issues.append(_issue(f"{path}/origin", "unknown channel origin", "metadata_value"))
    if "interpolation" in metadata and metadata["interpolation"] not in {"linear", "previous", "none"}:
        issues.append(_issue(f"{path}/interpolation", "unknown interpolation mode", "metadata_value"))


def _validate_channel_group(
    group: h5py.Group,
    path: str,
    issues: list[ValidationIssue],
    *,
    channel_id: str | None = None,
) -> tuple[Mapping[str, Any] | None, np.ndarray | None]:
    metadata_dataset = _check_dataset(group, "metadata_json", path, issues)
    metadata = _check_json_dataset(metadata_dataset, f"{path}/metadata_json", issues, object_name="channel metadata") if metadata_dataset is not None else None
    _metadata_fields(metadata, f"{path}/metadata_json", issues, channel_id=channel_id)
    values_dataset = _check_dataset(group, "values", path, issues)
    valid_dataset = _check_dataset(group, "valid", path, issues)
    values, _ = _check_values_and_mask(values_dataset, valid_dataset, f"{path}/values", f"{path}/valid", issues)
    return metadata, values


def _validate_track_group(group: h5py.Group, path: str, issues: list[ValidationIssue]) -> None:
    metadata_dataset = _check_dataset(group, "metadata_json", path, issues)
    _check_json_dataset(metadata_dataset, f"{path}/metadata_json", issues, object_name="track metadata") if metadata_dataset is not None else None
    distance_dataset = _check_dataset(group, "distance_m", path, issues)
    curvature_dataset = _check_dataset(group, "curvature_per_m", path, issues)
    distance = _check_numeric_vector(distance_dataset, f"{path}/distance_m", issues, name="track distance") if distance_dataset is not None else None
    curvature = _read_array(curvature_dataset, f"{path}/curvature_per_m", issues) if curvature_dataset is not None else None
    if curvature is not None:
        if curvature.ndim != 1:
            issues.append(_issue(f"{path}/curvature_per_m", "track curvature must be one-dimensional", "shape"))
        elif not np.issubdtype(curvature.dtype, np.number):
            issues.append(_issue(f"{path}/curvature_per_m", "track curvature must be numeric", "dtype"))
        else:
            curvature = np.asarray(curvature, dtype=np.float64)
            if not np.all(np.isfinite(curvature)):
                issues.append(_issue(f"{path}/curvature_per_m", "track curvature must contain only finite values", "finite"))
    if distance is not None and curvature is not None and distance.ndim == curvature.ndim == 1 and distance.shape != curvature.shape:
        issues.append(_issue(path, "track distance and curvature must have equal shapes", "shape"))
    x_dataset = group.get("x_m")
    y_dataset = group.get("y_m")
    if (x_dataset is None) != (y_dataset is None):
        issues.append(_issue(path, "track x_m and y_m must either both exist or both be absent", "geometry_pair"))
    geometry: dict[str, np.ndarray] = {}
    for name, dataset in (("x_m", x_dataset), ("y_m", y_dataset)):
        if dataset is None:
            continue
        if not isinstance(dataset, h5py.Dataset):
            issues.append(_issue(f"{path}/{name}", "track geometry must be a dataset", "dataset_type"))
            continue
        array = _read_array(dataset, f"{path}/{name}", issues)
        if array is None:
            continue
        geometry[name] = array
        if array.ndim != 1:
            issues.append(_issue(f"{path}/{name}", "track geometry must be one-dimensional", "shape"))
        elif not np.issubdtype(array.dtype, np.number):
            issues.append(_issue(f"{path}/{name}", "track geometry must be numeric", "dtype"))
        elif not np.all(np.isfinite(array)):
            issues.append(_issue(f"{path}/{name}", "track geometry must contain only finite values", "finite"))
        if distance is not None and array.ndim == 1 and array.shape != distance.shape:
            issues.append(_issue(f"{path}/{name}", "track geometry must match track distance shape", "shape"))


def _validate_axis_group(group: h5py.Group, path: str, issues: list[ValidationIssue]) -> None:
    time_dataset = _check_dataset(group, "time_s", path, issues)
    distance_dataset = _check_dataset(group, "distance_m", path, issues)
    time = _check_numeric_vector(time_dataset, f"{path}/time_s", issues, name="axis time") if time_dataset is not None else None
    distance = _check_numeric_vector(distance_dataset, f"{path}/distance_m", issues, name="axis distance") if distance_dataset is not None else None
    if time is not None and distance is not None and time.ndim == distance.ndim == 1 and time.shape != distance.shape:
        issues.append(_issue(path, "axis time and distance must have equal shapes", "shape"))


def _validate_lap_group(
    group: h5py.Group,
    path: str,
    issues: list[ValidationIssue],
    *,
    track_ids: set[str] | None = None,
) -> None:
    metadata_dataset = _check_dataset(group, "metadata_json", path, issues)
    lap_metadata: Mapping[str, Any] | None = None
    if metadata_dataset is not None:
        lap_metadata = _check_json_dataset(metadata_dataset, f"{path}/metadata_json", issues, object_name="lap metadata")
    if isinstance(lap_metadata, Mapping) and track_ids is not None:
        declared_track = lap_metadata.get("track_id") or lap_metadata.get("track")
        if declared_track is not None and str(declared_track) not in track_ids:
            issues.append(_issue(f"{path}/metadata_json/track_id", f"lap refers to unknown track {declared_track!r}", "unknown_track"))
    axes = group.get("axes")
    axis_lengths: dict[str, int] = {}
    if axes is None:
        issues.append(_issue(f"{path}/axes", "required group is missing", "missing_group"))
    elif not isinstance(axes, h5py.Group):
        issues.append(_issue(f"{path}/axes", "expected a group", "group_type"))
    else:
        axis_groups = [item for item in axes.values() if isinstance(item, h5py.Group)]
        if not axis_groups:
            issues.append(_issue(f"{path}/axes", "lap must contain at least one axis", "empty"))
        for axis_id, axis_group in axes.items():
            if isinstance(axis_group, h5py.Group):
                _validate_axis_group(axis_group, f"{path}/axes/{axis_id}", issues)
                time_dataset = axis_group.get("time_s")
                distance_dataset = axis_group.get("distance_m")
                if isinstance(time_dataset, h5py.Dataset) and isinstance(distance_dataset, h5py.Dataset):
                    try:
                        time = np.asarray(time_dataset[()])
                        distance = np.asarray(distance_dataset[()])
                    except Exception:
                        time = distance = np.asarray([])
                    if time.ndim == distance.ndim == 1 and time.shape == distance.shape:
                        axis_lengths[axis_id] = int(time.size)
            else:
                issues.append(_issue(f"{path}/axes/{axis_id}", "expected an axis group", "group_type"))
    channels = group.get("channels")
    if channels is None:
        issues.append(_issue(f"{path}/channels", "required group is missing", "missing_group"))
    elif not isinstance(channels, h5py.Group):
        issues.append(_issue(f"{path}/channels", "expected a group", "group_type"))
    else:
        for channel_id, channel_group in channels.items():
            if isinstance(channel_group, h5py.Group):
                channel_path = f"{path}/channels/{channel_id}"
                channel_metadata, values = _validate_channel_group(channel_group, channel_path, issues, channel_id=channel_id)
                if isinstance(channel_metadata, Mapping):
                    axis_id = channel_metadata.get("axis_id")
                    if axis_id is not None and not isinstance(axis_id, str):
                        issues.append(_issue(f"{channel_path}/metadata_json/axis_id", "channel axis_id must be a string", "metadata_value"))
                    elif axis_id is not None and axis_id not in axis_lengths:
                        issues.append(_issue(f"{channel_path}/metadata_json/axis_id", f"channel refers to unknown axis {axis_id!r}", "unknown_axis"))
                    elif axis_id in axis_lengths and values is not None and values.ndim == 1 and values.size != axis_lengths[axis_id]:
                        issues.append(_issue(f"{channel_path}/values", "channel values length does not match its axis", "shape"))
            else:
                issues.append(_issue(f"{path}/channels/{channel_id}", "expected a channel group", "group_type"))


def _validate_case_group(group: h5py.Group, path: str, issues: list[ValidationIssue], *, track_ids: set[str] | None = None) -> None:
    case_metadata: Mapping[str, Any] | None = None
    for name, object_name in (("metadata_json", "case metadata"), ("setup_json", "case setup")):
        dataset = _check_dataset(group, name, path, issues)
        if dataset is not None:
            payload = _check_json_dataset(dataset, f"{path}/{name}", issues, object_name=object_name)
            if name == "metadata_json" and isinstance(payload, Mapping):
                case_metadata = payload
    if case_metadata is not None and track_ids is not None:
        declared_track = case_metadata.get("track_id") or case_metadata.get("track")
        if declared_track is not None and str(declared_track) not in track_ids:
            issues.append(_issue(f"{path}/metadata_json/track_id", f"case refers to unknown track {declared_track!r}", "unknown_track"))
    laps = group.get("laps")
    if laps is None:
        issues.append(_issue(f"{path}/laps", "required group is missing", "missing_group"))
    elif not isinstance(laps, h5py.Group):
        issues.append(_issue(f"{path}/laps", "expected a group", "group_type"))
    else:
        for lap_id, lap_group in laps.items():
            if isinstance(lap_group, h5py.Group):
                _validate_lap_group(lap_group, f"{path}/laps/{lap_id}", issues, track_ids=track_ids)
            else:
                issues.append(_issue(f"{path}/laps/{lap_id}", "expected a lap group", "group_type"))
    envelope = group.get("envelope")
    if envelope is not None:
        if not isinstance(envelope, h5py.Group):
            issues.append(_issue(f"{path}/envelope", "expected a group", "group_type"))
        else:
            for channel_id, channel_group in envelope.items():
                if isinstance(channel_group, h5py.Group):
                    _validate_channel_group(channel_group, f"{path}/envelope/{channel_id}", issues, channel_id=channel_id)
                else:
                    issues.append(_issue(f"{path}/envelope/{channel_id}", "expected an envelope channel group", "group_type"))


def validate_file(path: str | Path) -> ValidationReport:
    """Validate an on-disk v1 HDF5 file without constructing result objects."""

    path = Path(path)
    issues: list[ValidationIssue] = []
    if not path.exists():
        return _report([_issue(str(path), "file does not exist", "missing_file")])
    if not path.is_file():
        return _report([_issue(str(path), "path is not a file", "file_type")])
    try:
        h5 = h5py.File(path, "r")
    except (OSError, ValueError) as error:
        return _report([_issue(str(path), f"could not open HDF5 file: {error}", "hdf5_open")])
    with h5:
        manifest_dataset = _check_dataset(h5, "manifest_json", "", issues)
        manifest = _check_json_dataset(manifest_dataset, "/manifest_json", issues, object_name="manifest") if manifest_dataset is not None else None
        if isinstance(manifest, Mapping):
            if manifest.get("format") != FORMAT_NAME:
                issues.append(_issue("/manifest_json/format", f"format must be {FORMAT_NAME!r}", "format"))
            version = manifest.get("schema_version")
            if version != f"{SCHEMA_VERSION[0]}.{SCHEMA_VERSION[1]}":
                issues.append(_issue("/manifest_json/schema_version", f"schema version must be {SCHEMA_VERSION[0]}.{SCHEMA_VERSION[1]!s}", "schema_version"))
            for key in ("run_uuid", "creation_utc", "source_type"):
                if key not in manifest:
                    issues.append(_issue(f"/manifest_json/{key}", "recommended manifest field is missing", "missing_manifest", "warning"))
        tracks = h5.get("tracks")
        track_ids: set[str] = set()
        if tracks is None:
            issues.append(_issue("/tracks", "required group is missing", "missing_group"))
        elif not isinstance(tracks, h5py.Group):
            issues.append(_issue("/tracks", "expected a group", "group_type"))
        else:
            for track_id, track_group in tracks.items():
                track_ids.add(track_id)
                if isinstance(track_group, h5py.Group):
                    _validate_track_group(track_group, f"/tracks/{track_id}", issues)
                else:
                    issues.append(_issue(f"/tracks/{track_id}", "expected a track group", "group_type"))
        cases = h5.get("cases")
        if cases is None:
            issues.append(_issue("/cases", "required group is missing", "missing_group"))
        elif not isinstance(cases, h5py.Group):
            issues.append(_issue("/cases", "expected a group", "group_type"))
        else:
            for case_id, case_group in cases.items():
                if isinstance(case_group, h5py.Group):
                    _validate_case_group(case_group, f"/cases/{case_id}", issues, track_ids=track_ids)
                else:
                    issues.append(_issue(f"/cases/{case_id}", "expected a case group", "group_type"))
    return _report(issues)


def _validate_axis(axis: Any, path: str, issues: list[ValidationIssue]) -> None:
    if not isinstance(axis, AxisData):
        issues.append(_issue(path, "expected AxisData", "type"))
        return
    time = np.asarray(axis.time_s)
    distance = np.asarray(axis.distance_m)
    if time.ndim != 1 or distance.ndim != 1 or time.shape != distance.shape or time.size == 0:
        issues.append(_issue(path, "axis time and distance must be equal-length non-empty vectors", "shape"))
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(distance)):
        issues.append(_issue(path, "axis coordinates must be finite", "finite"))
    if (time.size > 1 and np.any(np.diff(time) < 0)) or (distance.size > 1 and np.any(np.diff(distance) < 0)):
        issues.append(_issue(path, "axis coordinates must be monotonic non-decreasing", "monotonic"))


def _validate_channel(
    channel: Any,
    key: str,
    axes: Mapping[str, Any],
    path: str,
    issues: list[ValidationIssue],
    *,
    check_axis: bool = True,
) -> None:
    if not isinstance(channel, ChannelData):
        issues.append(_issue(path, "expected ChannelData", "type"))
        return
    if channel.metadata.id != key:
        issues.append(_issue(path, "channel map key must match metadata id", "metadata_id"))
    metadata = channel.metadata.as_dict()
    _metadata_fields(metadata, path, issues, channel_id=key)
    if check_axis and metadata.get("axis_id") not in axes:
        issues.append(_issue(f"{path}/axis_id", "channel refers to an unknown axis", "unknown_axis"))
    values = np.asarray(channel.values)
    valid = np.asarray(channel.valid)
    if values.ndim != 1 or valid.ndim != 1 or values.shape != valid.shape:
        issues.append(_issue(path, "channel values and valid mask must be equal-length vectors", "shape"))
        return
    if np.any(~valid & ~np.isnan(values)):
        issues.append(_issue(path, "invalid samples must be NaN", "nan_mask"))
    if np.any(valid & ~np.isfinite(values)):
        issues.append(_issue(path, "valid samples must be finite", "finite"))
    axis = axes.get(metadata.get("axis_id")) if check_axis else None
    if isinstance(axis, AxisData) and values.size != axis.time_s.size:
        issues.append(_issue(path, "channel length does not match its axis", "shape"))


def _validate_track(track: Any, key: str, path: str, issues: list[ValidationIssue]) -> None:
    if not isinstance(track, TrackData):
        issues.append(_issue(path, "expected TrackData", "type"))
        return
    if track.id != key:
        issues.append(_issue(path, "track map key must match track id", "id"))
    distance = np.asarray(track.distance_m)
    curvature = np.asarray(track.curvature_per_m)
    if distance.ndim != 1 or curvature.ndim != 1 or distance.size == 0 or distance.shape != curvature.shape:
        issues.append(_issue(path, "track distance and curvature must be equal-length non-empty vectors", "shape"))
    if not np.all(np.isfinite(distance)) or not np.all(np.isfinite(curvature)):
        issues.append(_issue(path, "track geometry arrays must be finite", "finite"))
    if distance.size > 1 and np.any(np.diff(distance) < 0):
        issues.append(_issue(path, "track distance must be monotonic non-decreasing", "monotonic"))
    if (track.x_m is None) != (track.y_m is None):
        issues.append(_issue(path, "track x_m and y_m must either both exist or both be absent", "geometry_pair"))
    for name in ("x_m", "y_m"):
        array = getattr(track, name)
        if array is not None:
            array = np.asarray(array)
            if array.ndim != 1 or array.shape != distance.shape:
                issues.append(_issue(f"{path}/{name}", "track geometry must match track distance shape", "shape"))
            elif not np.all(np.isfinite(array)):
                issues.append(_issue(f"{path}/{name}", "track geometry must be finite", "finite"))


def validate_lap(lap: LapData, track: TrackData | None = None) -> ValidationReport:
    """Return diagnostics for a :class:`LapData` and, optionally, its track."""

    issues: list[ValidationIssue] = []
    if not isinstance(lap, LapData):
        return _report([_issue("lap", "expected LapData", "type")])
    if not lap.axes:
        issues.append(_issue("lap/axes", "lap must contain at least one axis", "empty"))
    for axis_id, axis in lap.axes.items():
        _validate_axis(axis, f"lap/axes/{axis_id}", issues)
    for channel_id, channel in lap.channels.items():
        _validate_channel(channel, channel_id, lap.axes, f"lap/channels/{channel_id}", issues)
    if track is not None:
        _validate_track(track, lap.track_id, f"track/{lap.track_id}", issues)
    return _report(issues)


def validate_result_file(result: ResultFile) -> ValidationReport:
    """Return diagnostics for an already constructed :class:`ResultFile`."""

    issues: list[ValidationIssue] = []
    if not isinstance(result, ResultFile):
        return _report([_issue("result", "expected ResultFile", "type")])
    if result.manifest.get("format") != FORMAT_NAME:
        issues.append(_issue("manifest/format", f"format must be {FORMAT_NAME!r}", "format"))
    if result.manifest.get("schema_version") != f"{SCHEMA_VERSION[0]}.{SCHEMA_VERSION[1]}":
        issues.append(_issue("manifest/schema_version", "unsupported schema version", "schema_version"))
    for track_id, track in result.tracks.items():
        _validate_track(track, track_id, f"tracks/{track_id}", issues)
    for case_id, case in result.cases.items():
        if not isinstance(case, CaseData):
            issues.append(_issue(f"cases/{case_id}", "expected CaseData", "type"))
            continue
        for lap_id, lap in case.laps.items():
            if not isinstance(lap, LapData):
                issues.append(_issue(f"cases/{case_id}/laps/{lap_id}", "expected LapData", "type"))
                continue
            if lap.track_id not in result.tracks:
                issues.append(_issue(f"cases/{case_id}/laps/{lap_id}/track_id", f"lap refers to unknown track {lap.track_id!r}", "unknown_track"))
            report = validate_lap(lap, result.tracks.get(lap.track_id))
            issues.extend(
                ValidationIssue(f"cases/{case_id}/{issue.path}", issue.message, issue.code, issue.severity)
                for issue in report.issues
            )
        for channel_id, channel in case.envelope.items():
            _validate_channel(channel, channel_id, {}, f"cases/{case_id}/envelope/{channel_id}", issues, check_axis=False)
    return _report(issues)


# Friendly aliases used by integrations that refer to an HDF5 result rather
# than a generic file.
validate_h5 = validate_file
validate_hdf5 = validate_file
validate_results = validate_file


__all__ = [
    "ValidationIssue",
    "ValidationReport",
    "validate_file",
    "validate_h5",
    "validate_hdf5",
    "validate_lap",
    "validate_result_file",
    "validate_results",
]
