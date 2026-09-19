from __future__ import annotations

import numpy as np

from qss_telemetry.channel_catalog import channel_group, channel_role, default_channel_ids
from qss_telemetry.schema import ChannelData, ChannelMetadata


def _channel(channel_id: str, label: str | None = None) -> ChannelData:
    metadata = ChannelMetadata(
        id=channel_id,
        label=label or channel_id,
        unit="",
        axis_id="native",
        origin="simulation_output",
        interpolation="linear",
    )
    return ChannelData(metadata, np.array([0.0]), np.array([True]))


def test_catalog_classifies_scl_and_scd_as_lift_and_drag():
    assert channel_role("SCL") == "lift"
    assert channel_role("ClA") == "lift"
    assert channel_role("SCD") == "drag"
    assert channel_role("CdA") == "drag"
    assert channel_group("SCL") == "Aero"
    assert channel_group("SCD") == "Aero"


def test_default_channels_prioritize_driver_inputs_and_vehicle_constraints():
    channels = {
        channel_id: _channel(channel_id)
        for channel_id in (
            "misc_value",
            "speed_mps",
            "aSteer",
            "SCL",
            "SCD",
            "front_ride_height_in",
            "Fz_1",
            "throttle",
            "brake",
            "gear",
        )
    }

    selected = default_channel_ids(channels)

    assert selected[:2] == ["speed_mps", "aSteer"]
    assert {"SCL", "SCD", "front_ride_height_in", "Fz_1", "throttle", "brake", "gear"} <= set(selected)
    assert selected[-1] != "misc_value"
