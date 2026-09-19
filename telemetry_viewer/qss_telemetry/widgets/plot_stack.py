"""Linked, stacked pyqtgraph plots for one or two laps."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QLabel, QScrollArea, QVBoxLayout, QWidget

from ..channel_catalog import channel_role
from ..delta_mode import compute_delta_series
from ..schema import ChannelData, LapData
from .colors import ColorRegistry


def _axis_values(lap: LapData, channel: ChannelData, mode: str) -> np.ndarray:
    axis = lap.axes[channel.metadata.axis_id]
    return axis.time_s if mode == "time" else axis.distance_m


def _display_values(channel: ChannelData) -> np.ndarray:
    values = np.asarray(channel.values, dtype=float).copy()
    values[~np.asarray(channel.valid, dtype=bool)] = np.nan
    return values


class _ZoomOnlyViewBox(pg.ViewBox):
    """Keep wheel zoom while rejecting all mouse-drag panning."""

    pan_enabled = False

    def mouseDragEvent(self, ev, axis=None):  # noqa: N802 - pyqtgraph API
        del axis
        ev.ignore()


class PlotStack(QWidget):
    """Create one real PlotWidget per channel and link their x-ranges."""

    cursorValueChanged = Signal(float)

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
        self.datum_lap: LapData | None = None
        self.display_mode = "absolute"
        self.selected_channels: list[str] = []
        self.delta_errors: dict[str, str] = {}
        self._x_bounds: tuple[float, float] | None = None
        self._setting_cursor = False

    @property
    def plot_count(self) -> int:
        return len(self.plots)

    def set_data(
        self,
        reference_lap: LapData | None,
        selected_channels: Iterable[str],
        *,
        comparison_lap: LapData | None = None,
        datum_lap: LapData | None = None,
        display_mode: str = "absolute",
        axis_mode: str = "time",
        cursor_value: float | None = None,
        comparison_offset: float = 0.0,
        reference_offset: float = 0.0,
    ) -> None:
        self.reference_lap = reference_lap
        self.comparison_lap = comparison_lap
        self.datum_lap = datum_lap
        self.display_mode = "delta" if display_mode == "delta" else "absolute"
        self.selected_channels = list(dict.fromkeys(selected_channels))
        self.axis_mode = "distance" if axis_mode == "distance" else "time"
        self.comparison_offset = float(comparison_offset)
        self.reference_offset = float(reference_offset)
        self.delta_errors = {}
        self._x_bounds = None
        self._clear_plots()
        if reference_lap is None:
            self._set_empty("Select a lap to view telemetry")
            return
        if self.display_mode == "delta" and datum_lap is None:
            self._set_empty("Select a datum result to show delta waveforms")
            return
        plotted = 0
        for channel_id in self.selected_channels:
            ref_channel = reference_lap.channels.get(channel_id)
            comp_channel = comparison_lap.channels.get(channel_id) if comparison_lap is not None else None
            datum_channel = datum_lap.channels.get(channel_id) if datum_lap is not None else None
            if self.display_mode == "delta":
                if ref_channel is None or datum_channel is None:
                    continue
            elif ref_channel is None and comp_channel is None:
                continue
            channel_for_label = ref_channel or comp_channel or datum_channel
            if channel_for_label is None:
                continue
            plot = pg.PlotWidget(self.content, viewBox=_ZoomOnlyViewBox())
            plot.setBackground("#10131b")
            plot.showGrid(x=True, y=True, alpha=0.18)
            # The viewer owns horizontal navigation: wheel zoom and the
            # toolbar operate on every linked plot, while mouse dragging never
            # pans a single waveform away from the others.
            # x must remain enabled for pyqtgraph's wheelEvent to zoom.  The
            # custom view box above rejects the corresponding drag gestures.
            plot.getViewBox().setMouseEnabled(x=True, y=False)
            plot.setMinimumHeight(120)
            metadata = channel_for_label.metadata
            title = metadata.label or metadata.id
            role = channel_role(metadata.id, channel_for_label)
            if role in {"lift", "drag"}:
                title += f"  · {role.title()}"
            if self.display_mode == "delta":
                title = f"Delta {title} (active - datum)"
            elif comp_channel is None and comparison_lap is not None:
                title += "  · comparison missing"
            elif ref_channel is None and reference_lap is not None:
                title += "  · reference missing"
            plot.setTitle(title, color="#d9e2f2", size="10pt")
            plot.setLabel("left", metadata.label or metadata.id, units=metadata.unit or None)
            plot.setLabel("bottom", "Time", units="s") if self.axis_mode == "time" else plot.setLabel("bottom", "Distance", units="m")
            color = self.color_registry.color(channel_id)
            if self.display_mode == "delta":
                try:
                    series = compute_delta_series(reference_lap, datum_lap, channel_id, axis_mode=self.axis_mode)
                except (KeyError, ValueError) as error:
                    self.delta_errors[channel_id] = str(error)
                    continue
                x = series.x + float(reference_offset)
                curve = plot.plot(x, series.values, pen=pg.mkPen(color, width=2), name="Delta (active - datum)")
                curve.setClipToView(True)
                curve.setDownsampling(auto=True, method="peak")
                self.curves[(channel_id, "delta")] = curve
                self._include_x_values(x)
            else:
                if ref_channel is not None:
                    x = _axis_values(reference_lap, ref_channel, self.axis_mode) + float(reference_offset)
                    curve = plot.plot(x, _display_values(ref_channel), pen=pg.mkPen(color, width=2), name="Reference")
                    curve.setClipToView(True)
                    curve.setDownsampling(auto=True, method="peak")
                    self.curves[(channel_id, "reference")] = curve
                    self._include_x_values(x)
                if comp_channel is not None and comparison_lap is not None:
                    x = _axis_values(comparison_lap, comp_channel, self.axis_mode) + float(comparison_offset)
                    dashed = pg.mkPen(color, width=2, style=Qt.PenStyle.DashLine)
                    curve = plot.plot(x, _display_values(comp_channel), pen=dashed, name="Comparison")
                    curve.setClipToView(True)
                    curve.setDownsampling(auto=True, method="peak")
                    self.curves[(channel_id, "comparison")] = curve
                    self._include_x_values(x)
            if plotted == 0:
                plot.addLegend(offset=(8, 8))
            line = pg.InfiniteLine(angle=90, movable=True, pen=pg.mkPen("#f4c95d", width=1.4))
            line.sigPositionChanged.connect(self._cursor_line_moved)
            plot.addItem(line)
            self.plots.append(plot)
            self.cursor_lines.append(line)
            self.layout.insertWidget(self.layout.count() - 1, plot)
            if self.plots[:-1]:
                plot.setXLink(self.plots[0])
            plotted += 1
        if plotted == 0:
            self._set_empty("Selected channels are unavailable for this lap")
        else:
            if self._x_bounds is not None:
                for line in self.cursor_lines:
                    line.setBounds(self._x_bounds)
            self.fit_x_range()
            if cursor_value is not None:
                self.set_cursor_value(cursor_value)

    def set_cursor_value(self, value: float | None) -> None:
        if value is None:
            return
        self._setting_cursor = True
        try:
            value = float(value)
            if self._x_bounds is not None:
                value = float(np.clip(value, *self._x_bounds))
            for line in self.cursor_lines:
                line.setValue(value)
        finally:
            self._setting_cursor = False

    def fit_x_range(self) -> None:
        """Fit every linked plot to the complete reference/comparison domain."""

        if not self.plots or self._x_bounds is None:
            return
        self._set_linked_x_range(*self._x_bounds)

    def zoom_in(self, factor: float = 0.8) -> None:
        """Zoom in around the current shared center without changing y scales."""

        self._zoom_x(float(factor))

    def zoom_out(self, factor: float = 1.25) -> None:
        """Zoom out around the current shared center, capped at the fit range."""

        self._zoom_x(float(factor))

    def _zoom_x(self, factor: float) -> None:
        if not self.plots or factor <= 0.0:
            return
        current = self.plots[0].getViewBox().viewRange()[0]
        low, high = float(current[0]), float(current[1])
        center = 0.5 * (low + high)
        half_width = 0.5 * (high - low) * factor
        if half_width <= 0.0:
            return
        new_low = center - half_width
        new_high = center + half_width
        if self._x_bounds is not None:
            bound_low, bound_high = self._x_bounds
            new_low = max(bound_low, new_low)
            new_high = min(bound_high, new_high)
            if new_high <= new_low:
                new_low, new_high = bound_low, bound_high
        self._set_linked_x_range(new_low, new_high)

    def _set_linked_x_range(self, low: float, high: float) -> None:
        if high <= low:
            return
        # The x-link propagates this range, but setting all plots explicitly
        # also makes the invariant hold immediately in offscreen tests.
        for plot in self.plots:
            plot.setXRange(float(low), float(high), padding=0.0)

    def _include_x_values(self, values: np.ndarray) -> None:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return
        low, high = float(np.min(finite)), float(np.max(finite))
        if self._x_bounds is None:
            self._x_bounds = (low, high)
        else:
            self._x_bounds = (min(self._x_bounds[0], low), max(self._x_bounds[1], high))

    def _cursor_line_moved(self, line: pg.InfiniteLine) -> None:
        if self._setting_cursor:
            return
        value = float(line.value())
        if self._x_bounds is not None:
            value = float(np.clip(value, *self._x_bounds))
        self.set_cursor_value(value)
        self.cursorValueChanged.emit(value)

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
                datum_lap=self.datum_lap,
                display_mode=self.display_mode,
                axis_mode=mode,
                cursor_value=cursor_value,
                comparison_offset=self.comparison_offset,
                reference_offset=self.reference_offset,
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
