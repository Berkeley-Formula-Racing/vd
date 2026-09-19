"""Stable colours used by the telemetry viewer layers."""

from __future__ import annotations

import hashlib


# The order is deliberately fixed.  A channel keeps its colour when a user
# adds or removes a different channel, and the digest fallback keeps colours
# stable for channels added by a future exporter.
PALETTE = (
    "#4cc9f0",
    "#f72585",
    "#f8961e",
    "#90be6d",
    "#b5179e",
    "#577590",
    "#43aa8b",
    "#f94144",
    "#fee440",
    "#7209b7",
    "#4895ef",
    "#ff9f1c",
)


class ColorRegistry:
    """Map arbitrary channel IDs to deterministic display colours."""

    def __init__(self, palette: tuple[str, ...] = PALETTE) -> None:
        self.palette = palette or PALETTE
        self._colors: dict[str, str] = {}

    def color(self, key: str) -> str:
        if key not in self._colors:
            # A digest makes the mapping independent of selector order and
            # stable across sessions.  The palette is intentionally small and
            # high contrast; collisions are acceptable for unusually large
            # channel sets because line styles still identify the lap layer.
            digest = hashlib.sha1(key.encode("utf-8")).digest()
            color = self.palette[int.from_bytes(digest[:4], "big") % len(self.palette)]
            self._colors[key] = color
        return self._colors[key]

    def __getitem__(self, key: str) -> str:
        return self.color(key)

    def as_dict(self) -> dict[str, str]:
        return dict(self._colors)
