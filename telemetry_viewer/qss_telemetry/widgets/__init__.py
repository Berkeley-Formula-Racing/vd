"""Reusable Qt widgets for the QSS telemetry viewer."""

from .channel_selector import ChannelSelector, channel_group
from .colors import ColorRegistry
from .panels import DeltaDisplay, QualityPanel, SetupComparisonTable
from .playback import PlaybackController
from .plot_stack import PlotStack
from .track_map import TrackMapWidget, integrate_curvature

__all__ = [
    "ChannelSelector",
    "ColorRegistry",
    "DeltaDisplay",
    "PlaybackController",
    "PlotStack",
    "QualityPanel",
    "SetupComparisonTable",
    "TrackMapWidget",
    "channel_group",
    "integrate_curvature",
]
