"""Small deterministic v1 fixtures used by reader and UI tests."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

from .schema import AxisData, CaseData, ChannelData, ChannelMetadata, LapData, ResultFile, TrackData


def _json_bytes(value: object) -> np.ndarray:
    return np.frombuffer(json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8"), dtype=np.uint8)


def make_demo_result() -> ResultFile:
    """Return a four-sample lap with continuous and discrete channels."""
    distance_m = np.array([0.0, 10.0, 20.0, 30.0])
    time_s = np.array([0.0, 1.0, 2.0, 3.0])
    axis = AxisData("native", time_s, distance_m)
    channels = {
        "speed_mps": ChannelData(ChannelMetadata("speed_mps", "Vehicle speed", "m/s", "native", "simulation_output", "linear"), np.array([10.0, 10.0, 10.0, 10.0]), np.ones(4, dtype=bool)),
        "long_accel_mps2": ChannelData(ChannelMetadata("long_accel_mps2", "Longitudinal acceleration", "m/s^2", "native", "simulation_output", "linear"), np.zeros(4), np.ones(4, dtype=bool)),
        "lat_accel_mps2": ChannelData(ChannelMetadata("lat_accel_mps2", "Lateral acceleration", "m/s^2", "native", "simulation_output", "linear"), np.array([0.0, 2.0, 2.0, 0.0]), np.ones(4, dtype=bool)),
        "gear": ChannelData(ChannelMetadata("gear", "Gear", "1", "native", "simulation_output", "previous"), np.array([1.0, 1.0, 2.0, 2.0]), np.ones(4, dtype=bool)),
    }
    track = TrackData("demo_track", {"name": "Demo track", "geometry_source": "curvature"}, distance_m, np.array([0.0, 0.02, 0.02, 0.0]))
    lap = LapData("flying", {"role": "flying", "event": "autocross"}, {"native": axis}, channels, track.id)
    case = CaseData("baseline", {"label": "Baseline"}, {"mass_kg": 250.0}, {lap.id: lap})
    return ResultFile({"format": "qss-telemetry", "schema_version": "1.0", "run_uuid": "fixture-run", "source_type": "simulation"}, {track.id: track}, {case.id: case})


def write_demo_h5(path: str | Path) -> Path:
    """Write the demo result using the exact on-disk v1 layout."""
    result = make_demo_result()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(target, "w") as h5:
        h5.create_dataset("manifest_json", data=_json_bytes(result.manifest))
        track = result.tracks["demo_track"]
        track_group = h5.create_group("tracks/demo_track")
        track_group.create_dataset("metadata_json", data=_json_bytes(track.metadata))
        track_group.create_dataset("distance_m", data=track.distance_m)
        track_group.create_dataset("curvature_per_m", data=track.curvature_per_m)
        case = result.cases["baseline"]
        case_group = h5.create_group("cases/baseline")
        case_group.create_dataset("metadata_json", data=_json_bytes(case.metadata))
        case_group.create_dataset("setup_json", data=_json_bytes(case.setup))
        lap = case.laps["flying"]
        lap_group = case_group.create_group("laps/flying")
        lap_group.create_dataset("metadata_json", data=_json_bytes(lap.metadata))
        axis_group = lap_group.create_group("axes/native")
        axis_group.create_dataset("time_s", data=lap.axes["native"].time_s)
        axis_group.create_dataset("distance_m", data=lap.axes["native"].distance_m)
        for channel in lap.channels.values():
            group = lap_group.create_group(f"channels/{channel.metadata.id}")
            group.create_dataset("values", data=channel.values)
            group.create_dataset("valid", data=channel.valid.astype(np.uint8))
            group.create_dataset("metadata_json", data=_json_bytes(channel.metadata.as_dict()))
    return target
