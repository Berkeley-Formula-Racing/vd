import json

import h5py
import numpy as np
import pytest

from qss_telemetry.fixtures import write_demo_h5
from qss_telemetry.io import TelemetryFormatError, load_lap, load_results


def test_load_results_reads_the_complete_demo_tree_and_source_path(tmp_path):
    path = write_demo_h5(tmp_path / "demo.qss.h5")

    result = load_results(path)

    assert result.source_path == str(path)
    assert set(result.tracks) == {"demo_track"}
    assert set(result.cases) == {"baseline"}
    lap = load_lap(result, "baseline", "flying")
    assert set(lap.axes) == {"native"}
    assert set(lap.channels) == {"speed_mps", "long_accel_mps2", "lat_accel_mps2", "gear"}
    np.testing.assert_array_equal(lap.channels["gear"].valid, np.ones(4, dtype=bool))


def test_load_results_decodes_uint8_json_and_optional_track_geometry(tmp_path):
    path = write_demo_h5(tmp_path / "geometry.qss.h5")
    with h5py.File(path, "r+") as h5:
        track = h5["tracks/demo_track"]
        track.create_dataset("x_m", data=np.array([0.0, 10.0, 20.0, 30.0]))
        track.create_dataset("y_m", data=np.array([0.0, 1.0, 0.0, -1.0]))

    result = load_results(path)

    np.testing.assert_array_equal(result.tracks["demo_track"].x_m, [0.0, 10.0, 20.0, 30.0])
    np.testing.assert_array_equal(result.tracks["demo_track"].y_m, [0.0, 1.0, 0.0, -1.0])
    assert result.manifest["format"] == "qss-telemetry"


def test_load_results_reports_missing_required_dataset(tmp_path):
    path = write_demo_h5(tmp_path / "missing.qss.h5")
    with h5py.File(path, "r+") as h5:
        del h5["cases/baseline/laps/flying/channels/gear/valid"]

    with pytest.raises(TelemetryFormatError, match="valid") as error:
        load_results(path)

    assert not error.value.report.valid
    assert any(issue.path.endswith("/gear/valid") for issue in error.value.report.issues)


def test_load_results_rejects_incompatible_manifest_version(tmp_path):
    path = write_demo_h5(tmp_path / "version.qss.h5")
    with h5py.File(path, "r+") as h5:
        del h5["manifest_json"]
        payload = np.frombuffer(
            json.dumps({"format": "qss-telemetry", "schema_version": "2.0"}).encode("utf-8"),
            dtype=np.uint8,
        )
        h5.create_dataset("manifest_json", data=payload)

    with pytest.raises(TelemetryFormatError, match="schema_version"):
        load_results(path)
