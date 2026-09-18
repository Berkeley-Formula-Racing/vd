"""Track map widget with supplied-geometry and curvature fallback paths."""

from __future__ import annotations

from typing import Any

import numpy as np
import pyqtgraph as pg
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

from ..schema import TrackData


def integrate_curvature(track: TrackData) -> tuple[np.ndarray, np.ndarray]:
    """Integrate curvature into a planar schematic centerline.

    The first point is the origin and the initial heading points along +x.
    Trapezoidal heading integration keeps the result smooth when curvature is
    sampled sparsely and handles repeated distance samples without division.
    """

    distance = np.asarray(track.distance_m, dtype=float)
    curvature = np.asarray(track.curvature_per_m, dtype=float)
    x = np.zeros(distance.size, dtype=float)
    y = np.zeros(distance.size, dtype=float)
    if distance.size < 2:
        return x, y
    ds = np.diff(distance)
    heading = np.zeros(distance.size, dtype=float)
    heading[1:] = np.cumsum(0.5 * (curvature[:-1] + curvature[1:]) * ds)
    segment_heading = 0.5 * (heading[:-1] + heading[1:])
    x[1:] = np.cumsum(np.cos(segment_heading) * ds)
    y[1:] = np.cumsum(np.sin(segment_heading) * ds)
    return x, y


class TrackMapWidget(QWidget):
    """A real pyqtgraph map with one synchronized cursor marker."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.plot_widget = pg.PlotWidget(self)
        self.plot_widget.setBackground("#10131b")
        self.plot_widget.setAspectLocked(True)
        self.plot_widget.showGrid(x=True, y=True, alpha=0.15)
        self.plot_widget.setLabel("bottom", "x", units="m")
        self.plot_widget.setLabel("left", "y", units="m")
        self.quality_label = QLabel("No track loaded", self)
        self.quality_label.setObjectName("mapQualityLabel")
        self.quality_label.setStyleSheet("color: #f4c95d; font-weight: 600;")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.plot_widget, 1)
        layout.addWidget(self.quality_label)

        self.map_curve = pg.PlotDataItem(pen=pg.mkPen("#8aa4c8", width=2))
        self.marker_item = pg.ScatterPlotItem(size=13, brush=QColor("#ffca3a"), pen=pg.mkPen("#ffffff", width=1.2))
        self.plot_widget.addItem(self.map_curve)
        self.plot_widget.addItem(self.marker_item)
        self.track: TrackData | None = None
        self.track_distance_m = np.array([], dtype=float)
        self.x_m = np.array([], dtype=float)
        self.y_m = np.array([], dtype=float)
        self.marker_distance_m: float | None = None
        self.marker_position: tuple[float, float] | None = None
        self.has_supplied_geometry = False

    def set_track(self, track: TrackData | None) -> None:
        self.track = track
        if track is None:
            self.track_distance_m = np.array([], dtype=float)
            self.x_m = np.array([], dtype=float)
            self.y_m = np.array([], dtype=float)
            self.has_supplied_geometry = False
            self.map_curve.setData([], [])
            self.marker_item.setData([], [])
            self.marker_distance_m = None
            self.marker_position = None
            self.quality_label.setText("No track loaded")
            return
        self.track_distance_m = np.asarray(track.distance_m, dtype=float).copy()
        self.has_supplied_geometry = track.x_m is not None and track.y_m is not None
        if self.has_supplied_geometry:
            self.x_m = np.asarray(track.x_m, dtype=float).copy()
            self.y_m = np.asarray(track.y_m, dtype=float).copy()
            self.quality_label.setText("Supplied geometry")
        else:
            self.x_m, self.y_m = integrate_curvature(track)
            self.quality_label.setText("Schematic map")
        self.map_curve.setData(self.x_m, self.y_m)
        if self.x_m.size:
            self.plot_widget.setRange(
                xRange=(float(np.nanmin(self.x_m)), float(np.nanmax(self.x_m))),
                yRange=(float(np.nanmin(self.y_m)), float(np.nanmax(self.y_m))),
                padding=0.12,
            )
            self.set_marker_distance(float(self.track_distance_m[0]))
        else:
            self.marker_item.setData([], [])

    def set_marker_distance(self, distance_m: float | None) -> None:
        if distance_m is None or self.track_distance_m.size == 0:
            self.marker_distance_m = None
            self.marker_position = None
            self.marker_item.setData([], [])
            return
        bounded = float(np.clip(distance_m, self.track_distance_m[0], self.track_distance_m[-1]))
        x = float(np.interp(bounded, self.track_distance_m, self.x_m))
        y = float(np.interp(bounded, self.track_distance_m, self.y_m))
        self.marker_distance_m = bounded
        self.marker_position = (x, y)
        self.marker_item.setData([x], [y])

    def marker_data(self) -> dict[str, Any]:
        """Return a small serializable snapshot useful to integrations/tests."""

        return {
            "distance_m": self.marker_distance_m,
            "position": self.marker_position,
            "map_quality": self.quality_label.text(),
        }
