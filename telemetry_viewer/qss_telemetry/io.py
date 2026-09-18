"""Read version-1 QSS telemetry HDF5 files into the shared data contract."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np

from .schema import AxisData, CaseData, ChannelData, ChannelMetadata, LapData, ResultFile, SchemaError, TrackData
from .validation import ValidationReport, validate_file


class TelemetryFormatError(ValueError):
    """Raised when an HDF5 file does not satisfy the v1 interchange format."""

    def __init__(self, message: str, report: ValidationReport | None = None) -> None:
        super().__init__(message)
        self.report = report


def _json_payload(value: Any, path: str) -> bytes:
    """Convert an HDF5 scalar/vector payload into UTF-8 JSON bytes.

    The contract uses one-dimensional ``uint8`` datasets.  The bytes/string
    fallbacks make the reader tolerant of files produced by HDF5 wrappers
    that expose a byte vector as a scalar or a fixed-width string while the
    validator still reports those encodings as non-conforming.
    """

    array = np.asarray(value)
    if array.ndim == 1 and array.dtype == np.dtype(np.uint8):
        return bytes(array.tolist())
    if isinstance(value, bytes):
        return value
    if isinstance(value, str):
        return value.encode("utf-8")
    if array.ndim == 0:
        scalar = array.item()
        if isinstance(scalar, bytes):
            return scalar
        if isinstance(scalar, str):
            return scalar.encode("utf-8")
    if array.dtype.kind in {"S", "U"}:
        if array.dtype.kind == "S":
            return array.tobytes().rstrip(b"\x00")
        return str(array.item()).encode("utf-8")
    raise TelemetryFormatError(f"{path} must contain UTF-8 JSON bytes")


def _read_json(dataset: h5py.Dataset, path: str) -> Mapping[str, Any]:
    try:
        payload = _json_payload(dataset[()], path)
        value = json.loads(payload.decode("utf-8"))
    except TelemetryFormatError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError) as error:
        raise TelemetryFormatError(f"malformed JSON at {path}: {error}") from error
    if not isinstance(value, Mapping):
        raise TelemetryFormatError(f"JSON at {path} must contain an object")
    return dict(value)


def _required_dataset(group: h5py.Group, name: str, path: str) -> h5py.Dataset:
    item = group.get(name)
    if item is None:
        raise TelemetryFormatError(f"missing required dataset {path}/{name}")
    if not isinstance(item, h5py.Dataset):
        raise TelemetryFormatError(f"expected dataset at {path}/{name}")
    return item


def _array(group: h5py.Group, name: str, path: str, *, dtype: Any | None = None) -> np.ndarray:
    dataset = _required_dataset(group, name, path)
    try:
        array = np.asarray(dataset[()])
    except Exception as error:
        raise TelemetryFormatError(f"could not read dataset {path}/{name}: {error}") from error
    if dtype is not None:
        try:
            array = np.asarray(array, dtype=dtype)
        except (TypeError, ValueError, OverflowError) as error:
            raise TelemetryFormatError(f"dataset {path}/{name} is not {dtype}: {error}") from error
    return array


def _metadata(group: h5py.Group, path: str, name: str) -> Mapping[str, Any]:
    dataset = _required_dataset(group, name, path)
    return _read_json(dataset, f"{path}/{name}")


def _channel_metadata(payload: Mapping[str, Any], path: str) -> ChannelMetadata:
    required = ("id", "label", "unit", "axis_id", "origin", "interpolation")
    missing = [key for key in required if key not in payload]
    if missing:
        raise TelemetryFormatError(f"missing channel metadata at {path}: {', '.join(missing)}")
    try:
        return ChannelMetadata(
            id=str(payload["id"]),
            label=str(payload["label"]),
            unit=str(payload["unit"]),
            axis_id=str(payload["axis_id"]),
            origin=payload["origin"],
            interpolation=payload["interpolation"],
            description=str(payload.get("description", "")),
            coordinate_frame=str(payload.get("coordinate_frame", "vehicle")),
            sign_convention=str(payload.get("sign_convention", "")),
        )
    except (TypeError, ValueError, SchemaError) as error:
        raise TelemetryFormatError(f"invalid channel metadata at {path}: {error}") from error


def _load_channel(group: h5py.Group, channel_id: str, path: str) -> ChannelData:
    metadata = _channel_metadata(_metadata(group, path, "metadata_json"), f"{path}/metadata_json")
    values = _array(group, "values", path, dtype=np.float64)
    valid_raw = _array(group, "valid", path)
    if valid_raw.ndim != 1:
        raise TelemetryFormatError(f"channel validity mask at {path}/valid must be one-dimensional")
    if valid_raw.dtype != np.dtype(np.uint8) and valid_raw.dtype != np.dtype(bool):
        raise TelemetryFormatError(f"channel validity mask at {path}/valid must use uint8 values")
    if valid_raw.dtype == np.dtype(np.uint8) and np.any((valid_raw != 0) & (valid_raw != 1)):
        raise TelemetryFormatError(f"channel validity mask at {path}/valid must contain only 0 or 1")
    valid = np.asarray(valid_raw, dtype=bool)
    if channel_id != metadata.id:
        raise TelemetryFormatError(f"channel group {channel_id!r} does not match metadata id {metadata.id!r}")
    try:
        return ChannelData(metadata, values, valid)
    except SchemaError as error:
        raise TelemetryFormatError(f"invalid channel at {path}: {error}") from error


def _load_track(group: h5py.Group, track_id: str, path: str) -> TrackData:
    metadata = _metadata(group, path, "metadata_json")
    distance = _array(group, "distance_m", path, dtype=np.float64)
    curvature = _array(group, "curvature_per_m", path, dtype=np.float64)
    x_dataset = group.get("x_m")
    y_dataset = group.get("y_m")
    if (x_dataset is None) != (y_dataset is None):
        raise TelemetryFormatError(f"track geometry at {path} must include both x_m and y_m")
    x = np.asarray(x_dataset[()], dtype=np.float64) if isinstance(x_dataset, h5py.Dataset) else None
    y = np.asarray(y_dataset[()], dtype=np.float64) if isinstance(y_dataset, h5py.Dataset) else None
    try:
        return TrackData(track_id, metadata, distance, curvature, x, y)
    except SchemaError as error:
        raise TelemetryFormatError(f"invalid track at {path}: {error}") from error


def _load_axis(group: h5py.Group, axis_id: str, path: str) -> AxisData:
    time = _array(group, "time_s", path, dtype=np.float64)
    distance = _array(group, "distance_m", path, dtype=np.float64)
    try:
        return AxisData(axis_id, time, distance)
    except SchemaError as error:
        raise TelemetryFormatError(f"invalid axis at {path}: {error}") from error


def _load_lap(group: h5py.Group, lap_id: str, path: str, track_id: str) -> LapData:
    metadata = _metadata(group, path, "metadata_json")
    axes_group = group.get("axes")
    if not isinstance(axes_group, h5py.Group):
        raise TelemetryFormatError(f"missing required axes group at {path}/axes")
    axes: dict[str, AxisData] = {}
    for axis_id, axis_group in axes_group.items():
        if not isinstance(axis_group, h5py.Group):
            raise TelemetryFormatError(f"expected axis group at {path}/axes/{axis_id}")
        axes[axis_id] = _load_axis(axis_group, axis_id, f"{path}/axes/{axis_id}")
    channels_group = group.get("channels")
    if not isinstance(channels_group, h5py.Group):
        raise TelemetryFormatError(f"missing required channels group at {path}/channels")
    channels: dict[str, ChannelData] = {}
    for channel_id, channel_group in channels_group.items():
        if not isinstance(channel_group, h5py.Group):
            raise TelemetryFormatError(f"expected channel group at {path}/channels/{channel_id}")
        channels[channel_id] = _load_channel(channel_group, channel_id, f"{path}/channels/{channel_id}")
    try:
        return LapData(lap_id, metadata, axes, channels, track_id)
    except SchemaError as error:
        raise TelemetryFormatError(f"invalid lap at {path}: {error}") from error


def _manifest_case_track(manifest: Mapping[str, Any], case_id: str) -> Any | None:
    """Read a track association from common manifest inventory spellings."""

    for key in ("cases", "case_inventory", "case_inventory_json"):
        inventory = manifest.get(key)
        if isinstance(inventory, Mapping):
            entry = inventory.get(case_id)
            if isinstance(entry, Mapping):
                return entry.get("track_id") or entry.get("track")
        elif isinstance(inventory, (list, tuple)):
            for entry in inventory:
                if isinstance(entry, Mapping) and str(entry.get("id")) == case_id:
                    return entry.get("track_id") or entry.get("track")
    return None


def _load_case(
    group: h5py.Group,
    case_id: str,
    path: str,
    tracks: Mapping[str, TrackData],
    manifest: Mapping[str, Any],
) -> CaseData:
    metadata = _metadata(group, path, "metadata_json")
    setup = _metadata(group, path, "setup_json")
    case_track_id = metadata.get("track_id") or metadata.get("track") or group.attrs.get("track_id") or _manifest_case_track(manifest, case_id)
    laps_group = group.get("laps")
    if not isinstance(laps_group, h5py.Group):
        raise TelemetryFormatError(f"missing required laps group at {path}/laps")
    laps: dict[str, LapData] = {}
    for lap_id, lap_group in laps_group.items():
        if not isinstance(lap_group, h5py.Group):
            raise TelemetryFormatError(f"expected lap group at {path}/laps/{lap_id}")
        lap_metadata = _metadata(lap_group, f"{path}/laps/{lap_id}", "metadata_json")
        track_id = lap_metadata.get("track_id") or lap_metadata.get("track") or case_track_id
        if isinstance(track_id, (bytes, np.bytes_)):
            track_id = track_id.decode("utf-8")
        if track_id is None:
            # The fixture and the v1 tree carry the association as a case or
            # manifest inventory field in some exporters.  Fall back to the
            # sole track when unambiguous; otherwise fail with a useful path.
            if len(tracks) == 1:
                track_id = next(iter(tracks))
            else:
                raise TelemetryFormatError(f"lap {path}/laps/{lap_id} does not identify a track")
        track_id = str(track_id)
        if track_id not in tracks:
            raise TelemetryFormatError(f"lap {path}/laps/{lap_id} refers to unknown track {track_id!r}")
        laps[lap_id] = _load_lap(lap_group, lap_id, f"{path}/laps/{lap_id}", track_id)
    envelope_group = group.get("envelope")
    envelope: dict[str, ChannelData] = {}
    if envelope_group is not None:
        if not isinstance(envelope_group, h5py.Group):
            raise TelemetryFormatError(f"expected envelope group at {path}/envelope")
        for channel_id, channel_group in envelope_group.items():
            if not isinstance(channel_group, h5py.Group):
                raise TelemetryFormatError(f"expected envelope channel group at {path}/envelope/{channel_id}")
            envelope[channel_id] = _load_channel(channel_group, channel_id, f"{path}/envelope/{channel_id}")
    return CaseData(case_id, metadata, setup, laps, envelope)


def load_results(path: str | Path) -> ResultFile:
    """Load every track, case, lap, axis, and channel from a v1 HDF5 file.

    The returned object keeps ``source_path`` as the caller supplied path
    string.  Invalid files raise :class:`TelemetryFormatError`; the attached
    :attr:`TelemetryFormatError.report` contains all diagnostics found by
    :func:`validate_file`.
    """

    target = Path(path)
    report = validate_file(target)
    if not report.valid:
        raise TelemetryFormatError(f"invalid QSS telemetry file {target}: {report.summary}", report)
    try:
        with h5py.File(target, "r") as h5:
            manifest = dict(_read_json(_required_dataset(h5, "manifest_json", ""), "/manifest_json"))
            tracks_group = h5.get("tracks")
            tracks: dict[str, TrackData] = {}
            if isinstance(tracks_group, h5py.Group):
                for track_id, track_group in tracks_group.items():
                    if not isinstance(track_group, h5py.Group):
                        raise TelemetryFormatError(f"expected track group at /tracks/{track_id}")
                    tracks[track_id] = _load_track(track_group, track_id, f"/tracks/{track_id}")
            cases_group = h5.get("cases")
            cases: dict[str, CaseData] = {}
            if isinstance(cases_group, h5py.Group):
                for case_id, case_group in cases_group.items():
                    if not isinstance(case_group, h5py.Group):
                        raise TelemetryFormatError(f"expected case group at /cases/{case_id}")
                    cases[case_id] = _load_case(case_group, case_id, f"/cases/{case_id}", tracks, manifest)
    except TelemetryFormatError as error:
        if error.report is None:
            error.report = report
        raise
    except (OSError, ValueError, TypeError, SchemaError) as error:
        raise TelemetryFormatError(f"could not load QSS telemetry file {target}: {error}", report) from error
    return ResultFile(manifest, tracks, cases, source_path=str(path))


def load_lap(result_file: ResultFile, case_id: str, lap_id: str) -> LapData:
    """Return one lap from a loaded :class:`ResultFile` with a useful error."""

    if not isinstance(result_file, ResultFile):
        raise TypeError("result_file must be a ResultFile")
    return result_file.get_lap(case_id, lap_id)


__all__ = ["TelemetryFormatError", "load_lap", "load_results"]
