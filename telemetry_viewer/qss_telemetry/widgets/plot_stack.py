"""Linked, stacked pyqtgraph plots for one or two laps."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QScrollArea, QVBoxLayout, QWidget

from ..schema import ChannelData, LapData
from .colors import ColorRegistry


def _axis_values(lap: LapData, channel: ChannelData, mode: str) -> np.ndarray:
    axis = lap.axes[channel.metadata.axis_id]
    return axis.time_s if mode == "time" else axis.distance_m


def _display_values(channel: ChannelData) -> np.ndarray:
    values = np.asarray(channel.values, dtype=float).copy()
    values[~np.asarray(channel.valid, dtype=bool)] = np.nan
    return values


class PlotStack(QWidget):
    """Create one real PlotWidget per channel and link their x-ranges."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.scroll_area = QScrollArea(self)
        self.scroll_area.setWidgetResizable(True)
        self.content = QWidget()
        self.layout = QVBoxLayout(self.content)
        self.layout.setContentsMargins(2, 2, 2, 2)
        self.layout.setSpacing(4)
        self.layout.addStretch(1)
        self.scroll_area.setWidget(self.content)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(self.scroll_area)

        self.color_registry = ColorRegistry()
        self.plots: list[pg.PlotWidget] = []
        self.cursor_lines: list[pg.InfiniteLine] = []
        self.curves: dict[tuple[str, str], pg.PlotDataItem] = {}
        self._empty_label: QLabel | None = None
        self.axis_mode = "time"
        self.reference_lap: LapData | None = None
        self.comparison_lap: LapData | None = None
        self.selected_channels: list[str] = []

    @property
    def plot_count(self) -> int:
        return len(self.plots)

    def set_data(
        self,
        reference_lap: LapData | None,
        selected_channels: Iterable[str],
        *,
        comparison_lap: LapData | None = None,
        axis_mode: str = "time",
        cursor_value: float | None = None,
        comparison_offset: float = 0.0,
        reference_offset: float = 0.0,
    ) -> None:
        self.reference_lap = reference_lap
        self.comparison_lap = comparison_lap
        self.selected_channels = list(dict.fromkeys(selected_channels))
        self.axis_mode = "distance" if axis_mode == "distance" else "time"
        self._clear_plots()
        if reference_lap is None:
            self._set_empty("Select a lap to view telemetry")
            return
        plotted = 0
        for channel_id in self.selected_channels:
            ref_channel = reference_lap.channels.get(channel_id)
            comp_channel = comparison_lap.channels.get(channel_id) if comparison_lap is not None else None
            if ref_channel is None and comp_channel is None:
                continue
            channel_for_label = ref_channel or comp_channel
            if channel_for_label is None:
                continue
            plot = pg.PlotWidget(self.content)
            plot.setBackground("#10131b")
            plot.showGrid(x=True, y=True, alpha=0.18)
            plot.setMinimumHeight(120)
            metadata = channel_for_label.metadata
            title = metadata.label or metadata.id
            if comp_channel is None and comparison_lap is not None:
                title += "  · comparison missing"
            elif ref_channel is None and reference_lap is not None:
                title += "  · reference missing"
            plot.setTitle(title, color="#d9e2f2", size="10pt")
            plot.setLabel("left", metadata.label or metadata.id, units=metadata.unit or None)
            plot.setLabel("bottom", "Time", units="s") if self.axis_mode == "time" else plot.setLabel("bottom", "Distance", units="m")
            color = self.color_registry.color(channel_id)
            if ref_channel is not None:
                x = _axis_values(reference_lap, ref_channel, self.axis_mode) + float(reference_offset)
                curve = plot.plot(x, _display_values(ref_channel), pen=pg.mkPen(color, width=2), name="Reference")
                self.curves[(channel_id, "reference")] = curve
            if comp_channel is not None and comparison_lap is not None:
                x = _axis_values(comparison_lap, comp_channel, self.axis_mode) + float(comparison_offset)
                dashed = pg.mkPen(color, width=2, style=Qt.PenStyle.DashLine)
                curve = plot.plot(x, _display_values(comp_channel), pen=dashed, name="Comparison")
                self.curves[(channel_id, "comparison")] = curve
            if plotted == 0:
                plot.addLegend(offset=(8, 8))
            line = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen("#f4c95d", width=1.4))
            plot.addItem(line)
            self.plots.append(plot)
            self.cursor_lines.append(line)
            self.layout.insertWidget(self.layout.count() - 1, plot)
            if self.plots[:-1]:
                plot.setXLink(self.plots[0])
            plotted += 1
        if plotted == 0:
            self._set_empty("Selected channels are unavailable for this lap")
        elif cursor_value is not None:
            self.set_cursor_value(cursor_value)

    def set_cursor_value(self, value: float | None) -> None:
        if value is None:
            return
        for line in self.cursor_lines:
            line.setValue(float(value))

    def set_axis_mode(self, mode: str, cursor_value: float | None = None) -> None:
        mode = "distance" if mode == "distance" else "time"
        if mode == self.axis_mode and self.plots:
            self.set_cursor_value(cursor_value)
            return
        if self.reference_lap is not None:
            self.set_data(
                self.reference_lap,
                self.selected_channels,
                comparison_lap=self.comparison_lap,
                axis_mode=mode,
                cursor_value=cursor_value,
            )
        else:
            self.axis_mode = mode

    def _set_empty(self, message: str) -> None:
        self._empty_label = QLabel(message, self.content)
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #9aa9bd; padding: 24px;")
        self.layout.insertWidget(self.layout.count() - 1, self._empty_label)

    def _clear_plots(self) -> None:
        for plot in self.plots:
            self.layout.removeWidget(plot)
            plot.deleteLater()
        self.plots.clear()
        self.cursor_lines.clear()
        self.curves.clear()
        if self._empty_label is not None:
            self.layout.removeWidget(self._empty_label)
            self._empty_label.deleteLater()
            self._empty_label = None
