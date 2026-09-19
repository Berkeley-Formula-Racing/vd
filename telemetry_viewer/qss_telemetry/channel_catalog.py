"""Shared channel roles, groups, and default-view priorities.

Telemetry exporters use slightly different names for the same driver-facing
signals.  Keeping the matching rules here lets the model, selector, and any
future importer agree on what deserves to be visible first without changing
the on-disk channel IDs or units.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from .schema import ChannelData


@dataclass(frozen=True, slots=True)
class ChannelSpec:
    """A logical signal family understood by the viewer."""

    role: str
    group: str
    priority: int
    aliases: tuple[str, ...]
    contains: tuple[str, ...] = ()


_SPECS: tuple[ChannelSpec, ...] = (
    ChannelSpec("speed", "Vehicle dynamics", 10, ("speed_mps", "vcar", "vehicle_speed", "velocity"), ("speed",)),
    ChannelSpec("long_accel", "Vehicle dynamics", 11, ("long_accel_mps2", "glong", "longitudinal_acceleration"), ("long_accel",)),
    ChannelSpec("lat_accel", "Vehicle dynamics", 12, ("lat_accel_mps2", "glat", "lateral_acceleration"), ("lat_accel",)),
    ChannelSpec("yaw_rate", "Vehicle dynamics", 13, ("yaw_rate", "yaw_rate_radps"), ("yaw",)),
    ChannelSpec("steering", "Controls", 20, ("asteer", "steer_angle", "steering_angle", "steering_angle_rad", "steering"), ("steer",)),
    ChannelSpec("throttle", "Controls", 21, ("throttle", "throttle_demand", "accelerator"), ("throttle",)),
    ChannelSpec("brake", "Controls", 22, ("brake", "brakes", "brake_demand", "brake_pressure"), ("brake",)),
    ChannelSpec("gear", "Controls", 23, ("gear", "current_gear"), ("gear",)),
    ChannelSpec("control_demand", "Controls", 24, ("control_demand", "signed_control_demand", "control"), ("demand", "pedal")),
    # SCL/SCD are lift/drag channels.  ClA/CdA are accepted as exported
    # coefficient-area metrics; their source metadata remains authoritative.
    ChannelSpec("lift", "Aero", 30, ("scl", "cl", "cla", "lift", "lift_coefficient", "downforce", "downforce_n"), ("lift", "downforce", "cla")),
    ChannelSpec("drag", "Aero", 31, ("scd", "cd", "cda", "drag", "drag_coefficient", "drag_n"), ("drag", "cda")),
    ChannelSpec("ride_height", "Ride heights", 40, ("ride_height", "front_ride_height", "rear_ride_height"), ("ride_height", "rideheight")),
    ChannelSpec("tire_constraint", "Tire constraints", 50, ("fz", "fy", "fx", "alpha", "gamma", "kappa", "slip_angle", "slip_ratio"), ("tire", "tyre", "fz_", "fy_", "fx_", "alpha_", "gamma_", "kappa_")),
    ChannelSpec("lap_metric", "Lap metrics", 80, ("delta_time_s", "delta_s", "lap_delta_time_s", "curvature_per_m"), ("delta", "sector", "lap_time", "curvature")),
)


def _normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def _spec_for(channel_id: str) -> ChannelSpec | None:
    normalized = _normalize(channel_id)
    if not normalized:
        return None
    for spec in _SPECS:
        aliases = {_normalize(alias) for alias in spec.aliases}
        if normalized in aliases:
            return spec
        if any(normalized.startswith(f"{alias}_") for alias in aliases):
            return spec
        if any(token and token in normalized for token in spec.contains):
            return spec
    return None


def channel_spec(channel_id: str, channel: ChannelData | None = None) -> ChannelSpec | None:
    """Return the logical spec for a channel, if it is a recognized signal."""

    del channel  # Reserved for future metadata-aware matching.
    return _spec_for(channel_id)


def channel_role(channel_id: str, channel: ChannelData | None = None) -> str | None:
    """Return a stable logical role such as ``lift`` or ``drag``."""

    spec = channel_spec(channel_id, channel)
    return spec.role if spec is not None else None


def channel_group(channel_id: str, channel: ChannelData | None = None) -> str:
    """Return a friendly selector group, respecting exporter metadata first."""

    metadata = channel.metadata if channel is not None else None
    explicit = ""
    if metadata is not None:
        explicit = str(getattr(metadata, "group", "") or "")
        if not explicit and hasattr(metadata, "as_dict"):
            explicit = str(metadata.as_dict().get("group", "") or "")
    if explicit:
        return explicit.title()
    spec = channel_spec(channel_id, channel)
    return spec.group if spec is not None else "Other"


def channel_priority(channel_id: str, channel: ChannelData | None = None) -> int | None:
    """Return the default-view priority, or ``None`` for unclassified data."""

    spec = channel_spec(channel_id, channel)
    return spec.priority if spec is not None else None


def default_channel_ids(channels: Mapping[str, ChannelData] | Iterable[ChannelData]) -> list[str]:
    """Return recognized channels in a useful driver-analysis order.

    Unknown channels remain available in the selector but are intentionally
    not added to the initial plot stack.  This keeps a large result responsive
    while still exposing every exporter-provided signal on demand.
    """

    values = list(channels.values()) if isinstance(channels, Mapping) else list(channels)
    ranked: list[tuple[int, str, str]] = []
    for channel in values:
        channel_id = str(channel.metadata.id)
        priority = channel_priority(channel_id, channel)
        if priority is None:
            continue
        ranked.append((priority, channel.metadata.label.casefold(), channel_id))
    ranked.sort()
    return list(dict.fromkeys(item[2] for item in ranked))


__all__ = [
    "ChannelSpec",
    "channel_group",
    "channel_priority",
    "channel_role",
    "channel_spec",
    "default_channel_ids",
]
