from qss_telemetry.fixtures import make_demo_result, write_demo_h5


def test_demo_result_obeys_shared_data_contract(tmp_path):
    result = make_demo_result()
    lap = result.get_lap("baseline", "flying")

    assert lap.axes["native"].distance_m[-1] == 30.0
    assert lap.channels["gear"].metadata.interpolation == "previous"
    assert write_demo_h5(tmp_path / "demo.qss.h5").is_file()
