import json

import h5py
import numpy as np

from qss_telemetry.fixtures import make_demo_result, write_demo_h5
from qss_telemetry.io import load_results
from qss_telemetry.validation import validate_file, validate_lap


def test_validate_lap_accepts_demo_lap():
    report = validate_lap(make_demo_result().get_lap("baseline", "flying"))

    assert report.valid
    assert report.issues == ()


def test_validate_file_reports_nonmonotonic_axis_and_invalid_mask_values(tmp_path):
    path = write_demo_h5(tmp_path / "corrupt.qss.h5")
    with h5py.File(path, "r+") as h5:
        h5["cases/baseline/laps/flying/axes/native/time_s"][...] = [0.0, 2.0, 1.0, 3.0]
        h5["cases/baseline/laps/flying/channels/speed_mps/values"][...] = [10.0, np.nan, 10.0, 10.0]
        h5["cases/baseline/laps/flying/channels/speed_mps/valid"][...] = [1, 1, 1, 1]

    report = validate_file(path)

    assert not report.valid
    assert any("monotonic" in issue.message for issue in report.issues)
    assert any("finite" in issue.message or "NaN" in issue.message for issue in report.issues)


def test_validate_file_reports_malformed_geometry_pair(tmp_path):
    path = write_demo_h5(tmp_path / "geometry-corrupt.qss.h5")
    with h5py.File(path, "r+") as h5:
        h5["tracks/demo_track"].create_dataset("x_m", data=np.array([0.0, 1.0, 2.0]))

    report = validate_file(path)

    assert not report.valid
    assert any("x_m" in issue.path or "geometry" in issue.message for issue in report.issues)


def test_validate_file_reports_unknown_channel_axis_and_length_mismatch(tmp_path):
    path = write_demo_h5(tmp_path / "channel-shape-corrupt.qss.h5")
    with h5py.File(path, "r+") as h5:
        metadata_path = "cases/baseline/laps/flying/channels/speed_mps/metadata_json"
        del h5[metadata_path]
        metadata = {
            "id": "speed_mps",
            "label": "Vehicle speed",
            "unit": "m/s",
            "axis_id": "missing-axis",
            "origin": "simulation_output",
            "interpolation": "linear",
            "description": "",
            "coordinate_frame": "vehicle",
            "sign_convention": "",
        }
        h5.create_dataset(metadata_path, data=np.frombuffer(json.dumps(metadata).encode("utf-8"), dtype=np.uint8))
        del h5["cases/baseline/laps/flying/channels/long_accel_mps2/values"]
        h5["cases/baseline/laps/flying/channels/long_accel_mps2"].create_dataset("values", data=np.zeros(3))

    report = validate_file(path)

    assert not report.valid
    assert any(issue.code == "unknown_axis" for issue in report.issues)
    assert any(issue.code == "shape" and "long_accel" in issue.path for issue in report.issues)
