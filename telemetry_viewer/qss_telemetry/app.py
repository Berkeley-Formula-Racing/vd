"""Local desktop viewer for versioned QSS telemetry results.

The application consumes :class:`~qss_telemetry.schema.ResultFile` objects
directly.  The HDF5 reader is imported only by :func:`main` and the Open-file
action, which keeps the UI usable with the deterministic fixture and gives a
small adapter boundary for future result sources.
"""

from __future__ import annotations

import csv
import math
import sys
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from .schema import AxisData, CaseData, ChannelData, LapData, ResultFile, TrackData


ResultInput = ResultFile | Mapping[str, ResultFile] | Iterable[ResultFile]
DeltaTimeProvider = Callable[[LapData, LapData], float | Mapping[str, Any] | None]


def _result_map(value: ResultInput | None) -> dict[str, ResultFile]:
    if value is None:
        return {}
    if isinstance(value, ResultFile):
        name = Path(value.source_path).stem if value.source_path else "Result"
        return {name or "Result": value}
    if isinstance(value, Mapping):
        return {str(key): item for key, item in value.items() if isinstance(item, ResultFile)}
    return {f"Result {index + 1}": item for index, item in enumerate(value) if isinstance(item, ResultFile)}


def _first_case_lap(result: ResultFile | None) -> tuple[str | None, str | None]:
    if result is None or not result.cases:
        return None, None
    case_id, case = next(iter(result.cases.items()))
    lap_id = next(iter(case.laps), None)
    return case_id, lap_id


def _preferred_axis(lap: LapData | None) -> AxisData | None:
    if lap is None or not lap.axes:
        return None
    if "native" in lap.axes:
        return lap.axes["native"]
    return next(iter(lap.axes.values()))


def _axis_coordinates(lap: LapData | None) -> tuple[np.ndarray, np.ndarray]:
    axis = _preferred_axis(lap)
    if axis is None:
        return np.array([], dtype=float), np.array([], dtype=float)
    return axis.time_s, axis.distance_m


def integrate_curvature(track: TrackData) -> tuple[np.ndarray, np.ndarray]:
    """Integrate a curvature-only track into a clearly schematic centerline."""

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


def track_geometry(track: TrackData | None) -> tuple[np.ndarray, np.ndarray, str]:
    """Return x/y geometry and the quality label shown in the UI."""

    if track is None:
        return np.array([], dtype=float), np.array([], dtype=float), "No track"
    if track.x_m is not None and track.y_m is not None:
        return np.asarray(track.x_m, dtype=float), np.asarray(track.y_m, dtype=float), "Supplied geometry"
    x, y = integrate_curvature(track)
    return x, y, "Schematic map"


def _delta_time_result(
    reference_lap: LapData | None,
    comparison_lap: LapData | None,
    provider: DeltaTimeProvider | None = None,
    alignment_settings: Mapping[str, Any] | None = None,
) -> tuple[float | None, str | None]:
    """Return ``(delta, error)`` using the core overlap semantics when present."""

    if reference_lap is None or comparison_lap is None:
        return None, None
    if provider is not None:
        try:
            result = provider(reference_lap, comparison_lap)
            if isinstance(result, Mapping):
                for key in ("delta_time_s", "delta_seconds", "delta_time", "dt_s"):
                    if key in result and result[key] is not None:
                        return float(result[key]), None
                return None, "delta provider returned no delta-time value"
            if result is not None:
                return float(result), None
            return None, "delta provider returned no value"
        except Exception as error:
            return None, str(error)

    try:
        # The core package exposes alignment as ``align_laps`` followed by
        # ``calculate_delta_time(AlignedPair)``.  Import it lazily so the
        # fixture driven UI remains independent while an integrated build gets
        # the exact shared overlap/gap semantics automatically.
        from .alignment import align_laps, calculate_delta_time

        aligned = align_laps(reference_lap, comparison_lap, alignment_settings)
        delta_channel = calculate_delta_time(aligned)
        valid = np.asarray(delta_channel.valid, dtype=bool)
        values = np.asarray(delta_channel.values, dtype=float)
        if np.any(valid):
            return float(values[np.flatnonzero(valid)[-1]]), None
        return None, "laps have no common overlap"
    except ImportError:
        # A standalone fixture build may omit the optional core module.  Its
        # simple duration fallback is useful for development in that mode.
        pass
    except Exception as error:
        # Once the core is importable, alignment failures are meaningful: a
        # mismatched track or empty overlap must remain visible to the user.
        return None, str(error)

    reference_time, _ = _axis_coordinates(reference_lap)
    comparison_time, _ = _axis_coordinates(comparison_lap)
    if reference_time.size == 0 or comparison_time.size == 0:
        return None, "selected laps have no time samples"
    return float((comparison_time[-1] - comparison_time[0]) - (reference_time[-1] - reference_time[0])), None


def delta_time_seconds(
    reference_lap: LapData | None,
    comparison_lap: LapData | None,
    provider: DeltaTimeProvider | None = None,
    alignment_settings: Mapping[str, Any] | None = None,
) -> float | None:
    """Return comparison minus reference duration, positive when slower.

    ``provider`` is the injection point for the core alignment implementation.
    It may return a number or a mapping containing one of the usual delta-time
    keys; when it is absent, the displayed value is derived from each lap's
    native axis duration.
    """

    value, _ = _delta_time_result(reference_lap, comparison_lap, provider, alignment_settings)
    return value


def _channel_group(channel_id: str, channel: ChannelData | None = None) -> str:
    metadata = channel.metadata if channel is not None else None
    value = channel_id.lower()
    label = metadata.label.lower() if metadata is not None else ""
    text = f"{value} {label}"
    if any(token in text for token in ("speed", "velocity", "accel", "yaw", "slip", "kinematic")):
        return "Vehicle dynamics"
    if any(token in text for token in ("steer", "throttle", "brake", "control", "demand", "pedal", "torque")):
        return "Controls"
    if any(token in text for token in ("wheel", "tire", "tyre", "suspension", "damper", "load")):
        return "Wheels and suspension"
    if any(token in text for token in ("delta", "lap", "sector", "time")):
        return "Lap metrics"
    return "Other"


class ViewerModel:
    """Qt-independent state and adapters for the desktop viewer."""

    def __init__(self, results: ResultInput | None = None, *, delta_time_provider: DeltaTimeProvider | None = None) -> None:
        self.results = _result_map(results)
        self.result_name = next(iter(self.results), None)
        self.case_id, self.lap_id = _first_case_lap(self.current_result)
        self.reference_lap_id = self.lap_id
        self.comparison_lap_id: str | None = None
        self.axis_mode = "time"
        self.distance_offset_m = 0.0
        self.time_offset_s = 0.0
        self.cursor_index = 0
        self.selected_channel_ids: list[str] = []
        self._missing_channels: list[str] = []
        self._comparison_error: str | None = None
        self.delta_time_provider = delta_time_provider
        self._select_default_channels()

    @property
    def current_result(self) -> ResultFile | None:
        return self.results.get(self.result_name) if self.result_name is not None else None

    @property
    def current_case(self) -> CaseData | None:
        result = self.current_result
        return result.cases.get(self.case_id) if result is not None and self.case_id is not None else None

    @property
    def current_lap(self) -> LapData | None:
        case = self.current_case
        return case.laps.get(self.lap_id) if case is not None and self.lap_id is not None else None

    @property
    def reference_lap(self) -> LapData | None:
        case = self.current_case
        if case is None or self.reference_lap_id is None:
            return None
        return case.laps.get(self.reference_lap_id)

    @property
    def comparison_lap(self) -> LapData | None:
        case = self.current_case
        if case is None or self.comparison_lap_id is None:
            return None
        return case.laps.get(self.comparison_lap_id)

    @property
    def available_channels(self) -> dict[str, ChannelData]:
        channels: dict[str, ChannelData] = {}
        for lap in (self.reference_lap, self.comparison_lap):
            if lap is not None:
                channels.update(lap.channels)
        return channels

    @property
    def missing_channels(self) -> list[str]:
        return list(self._missing_channels)

    @property
    def track(self) -> TrackData | None:
        result = self.current_result
        lap = self.reference_lap or self.current_lap
        if result is None or lap is None:
            return None
        track = result.tracks.get(lap.track_id)
        return track if track is not None else next(iter(result.tracks.values()), None)

    @property
    def cursor_time_s(self) -> float | None:
        axis = _preferred_axis(self.reference_lap or self.current_lap)
        if axis is None:
            return None
        index = min(max(0, self.cursor_index), axis.time_s.size - 1)
        return float(axis.time_s[index])

    @property
    def cursor_distance_m(self) -> float | None:
        axis = _preferred_axis(self.reference_lap or self.current_lap)
        if axis is None:
            return None
        index = min(max(0, self.cursor_index), axis.distance_m.size - 1)
        return float(axis.distance_m[index])

    @property
    def cursor_value(self) -> float | None:
        return self.cursor_distance_m if self.axis_mode == "distance" else self.cursor_time_s

    @property
    def delta_seconds(self) -> float | None:
        value, error = _delta_time_result(
            self.reference_lap,
            self.comparison_lap,
            self.delta_time_provider,
            {
                # Delta time is a physical-distance comparison even when the
                # plot x-axis is switched to time.
                "mode": "distance",
                "distance_offset_m": self.distance_offset_m,
                "time_offset_s": self.time_offset_s,
            },
        )
        self._comparison_error = error
        return value

    @property
    def comparison_error(self) -> str | None:
        # Evaluate once so status consumers always see the latest offsets.
        self.delta_seconds
        return self._comparison_error

    @property
    def comparison_status(self) -> str:
        if self.comparison_lap is None:
            return "No comparison selected"
        if self.delta_seconds is None:
            return "Comparison unavailable: " + (self._comparison_error or "laps have no common overlap")
        return "Comparison aligned"

    @property
    def setup_comparison(self) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
        reference = self.reference_lap
        comparison = self.comparison_lap
        reference_setup = self.current_case.setup if self.current_case is not None else {}
        # A future exporter may persist lap-level setup overrides.
        if reference is not None and isinstance(reference.metadata.get("setup"), Mapping):
            reference_setup = reference.metadata["setup"]
        comparison_setup = reference_setup
        if comparison is not None and isinstance(comparison.metadata.get("setup"), Mapping):
            comparison_setup = comparison.metadata["setup"]
        return reference_setup, comparison_setup

    @property
    def quality_summary(self) -> dict[str, Any]:
        lap = self.reference_lap or self.current_lap
        track = self.track
        total = sum(channel.valid.size for channel in lap.channels.values()) if lap is not None else 0
        invalid = sum(int((~channel.valid).sum()) for channel in lap.channels.values()) if lap is not None else 0
        metadata = dict(lap.metadata) if lap is not None else {}
        reconstruction = metadata.get("reconstruction") or metadata.get("qss_reconstruction") or metadata.get("reconstruction_status")
        return {
            "map_quality": "Supplied geometry" if track is not None and track.x_m is not None else ("Schematic map" if track is not None else "No track"),
            "reconstruction": reconstruction if reconstruction is not None else "not reported",
            "missing_channels": self.missing_channels,
            "comparison_status": self.comparison_status,
            "valid_samples": total - invalid,
            "invalid_samples": invalid,
            "sample_count": total,
        }

    def select_result(self, result_name: str) -> None:
        if result_name not in self.results:
            raise KeyError(f"unknown result: {result_name}")
        self.result_name = result_name
        self.case_id, self.lap_id = _first_case_lap(self.current_result)
        self.reference_lap_id = self.lap_id
        self.comparison_lap_id = None
        self._comparison_error = None
        self.cursor_index = 0
        self._select_default_channels()

    def select_lap(self, case_id: str, lap_id: str) -> None:
        result = self.current_result
        if result is None or case_id not in result.cases or lap_id not in result.cases[case_id].laps:
            raise KeyError(f"unknown case/lap: {case_id}/{lap_id}")
        self.case_id = case_id
        self.lap_id = lap_id
        self.reference_lap_id = lap_id
        self.comparison_lap_id = None
        self._comparison_error = None
        self.cursor_index = 0
        self._select_default_channels()

    def set_comparison_laps(self, reference_lap_id: str, comparison_lap_id: str | None, case_id: str | None = None) -> None:
        if case_id is not None and case_id != self.case_id:
            self.select_lap(case_id, reference_lap_id)
        case = self.current_case
        if case is None or reference_lap_id not in case.laps:
            raise KeyError(f"unknown reference lap: {reference_lap_id}")
        if comparison_lap_id is not None and comparison_lap_id not in case.laps:
            raise KeyError(f"unknown comparison lap: {comparison_lap_id}")
        self.reference_lap_id = reference_lap_id
        self.lap_id = reference_lap_id
        self.comparison_lap_id = comparison_lap_id
        self._comparison_error = None
        self.cursor_index = min(self.cursor_index, max(0, len(_preferred_axis(self.reference_lap).time_s) - 1)) if _preferred_axis(self.reference_lap) is not None else 0
        self._select_default_channels(keep_existing=True)

    def set_selected_channels(self, channel_ids: Iterable[str]) -> list[str]:
        requested = list(dict.fromkeys(str(item) for item in channel_ids))
        available = self.available_channels
        self._missing_channels = [item for item in requested if item not in available]
        self.selected_channel_ids = [item for item in requested if item in available]
        return self.missing_channels

    def set_axis_mode(self, mode: str) -> None:
        self.axis_mode = "distance" if mode == "distance" else "time"

    def set_alignment_offsets(self, *, distance_offset_m: float | None = None, time_offset_s: float | None = None) -> None:
        if distance_offset_m is not None:
            self.distance_offset_m = float(distance_offset_m)
        if time_offset_s is not None:
            self.time_offset_s = float(time_offset_s)
        self._comparison_error = None

    def set_cursor_index(self, index: int) -> int:
        axis = _preferred_axis(self.reference_lap or self.current_lap)
        maximum = max(0, axis.time_s.size - 1) if axis is not None else 0
        self.cursor_index = min(max(0, int(index)), maximum)
        return self.cursor_index

    def set_cursor_value(self, value: float) -> int:
        axis = _preferred_axis(self.reference_lap or self.current_lap)
        if axis is None or axis.time_s.size == 0:
            self.cursor_index = 0
            return self.cursor_index
        values = axis.distance_m if self.axis_mode == "distance" else axis.time_s
        self.cursor_index = int(np.argmin(np.abs(values - float(value))))
        return self.cursor_index

    def grouped_channels(self) -> dict[str, list[ChannelData]]:
        groups: dict[str, list[ChannelData]] = {}
        for channel_id, channel in self.available_channels.items():
            groups.setdefault(_channel_group(channel_id, channel), []).append(channel)
        for values in groups.values():
            values.sort(key=lambda item: (item.metadata.label.casefold(), item.metadata.id))
        return dict(sorted(groups.items(), key=lambda pair: pair[0]))

    def _select_default_channels(self, *, keep_existing: bool = False) -> None:
        available = self.available_channels
        if keep_existing and self.selected_channel_ids:
            self.set_selected_channels(self.selected_channel_ids)
            return
        preferred: list[str] = []
        exact_order = (
            "speed_mps",
            "long_accel_mps2",
            "lat_accel_mps2",
            "control",
            "control_demand",
            "steering",
            "steering_angle_rad",
            "delta_time_s",
            "delta_s",
            "lap_delta_time_s",
            # The fixture has gear in place of control/delta, keeping the
            # default view useful while preserving the required core channels.
            "gear",
        )
        for channel_id in exact_order:
            if channel_id in available and channel_id not in preferred:
                preferred.append(channel_id)
        for channel_id, channel in available.items():
            text = f"{channel_id} {channel.metadata.label}".lower()
            if any(token in text for token in ("control", "demand", "delta")) and channel_id not in preferred:
                preferred.append(channel_id)
        if not preferred:
            preferred = list(available)[:4]
        self.set_selected_channels(preferred)


try:  # Keep ViewerModel importable for non-Qt tooling and headless validation.
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import (
        QApplication,
        QComboBox,
        QDoubleSpinBox,
        QFileDialog,
        QHBoxLayout,
        QLabel,
        QMainWindow,
        QPushButton,
        QSlider,
        QSplitter,
        QTabWidget,
        QTreeWidget,
        QTreeWidgetItem,
        QVBoxLayout,
        QWidget,
    )
    import pyqtgraph as pg

    from .widgets import ChannelSelector, DeltaDisplay, PlaybackController, PlotStack, QualityPanel, SetupComparisonTable, TrackMapWidget

    _QT_AVAILABLE = True
    _QT_IMPORT_ERROR: Exception | None = None
except Exception as error:  # pragma: no cover - exercised only without UI dependencies
    _QT_AVAILABLE = False
    _QT_IMPORT_ERROR = error


if _QT_AVAILABLE:

    class TelemetryViewerWindow(QMainWindow):
        """Interactive QSS result viewer with an injectable model."""

        def __init__(
            self,
            result: ResultInput | None = None,
            *,
            results: ResultInput | None = None,
            parent: QWidget | None = None,
            delta_time_provider: DeltaTimeProvider | None = None,
        ) -> None:
            super().__init__(parent)
            self.setWindowTitle("QSS Telemetry Viewer")
            self.resize(1440, 900)
            self.model = ViewerModel(result if result is not None else results, delta_time_provider=delta_time_provider)
            self._build_ui()
            self._populate_browser()
            self._refresh_view(reset_playback=True)

        @property
        def cursor_index(self) -> int:
            return self.model.cursor_index

        @property
        def delta_seconds(self) -> float | None:
            return self.model.delta_seconds

        @property
        def missing_channels(self) -> list[str]:
            return self.model.missing_channels

        @property
        def reference_lap(self) -> LapData | None:
            return self.model.reference_lap

        @property
        def comparison_lap(self) -> LapData | None:
            return self.model.comparison_lap

        def _build_ui(self) -> None:
            central = QWidget(self)
            self.setCentralWidget(central)
            root = QVBoxLayout(central)
            root.setContentsMargins(8, 8, 8, 8)
            root.setSpacing(6)

            toolbar = QHBoxLayout()
            self.open_button = QPushButton("Open result…", self)
            self.export_csv_button = QPushButton("Export CSV…", self)
            self.export_screenshot_button = QPushButton("Screenshot…", self)
            self.axis_mode_combo = QComboBox(self)
            self.axis_mode_combo.addItem("Time", "time")
            self.axis_mode_combo.addItem("Distance", "distance")
            self.play_button = QPushButton("Play", self)
            self.speed_combo = QComboBox(self)
            for speed in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
                self.speed_combo.addItem(f"{speed:g}×", speed)
            self.speed_combo.setCurrentIndex(2)
            self.cursor_time_label = QLabel("t=— · s=—", self)
            self.delta_label = DeltaDisplay(self)
            toolbar.addWidget(self.open_button)
            toolbar.addWidget(self.export_csv_button)
            toolbar.addWidget(self.export_screenshot_button)
            toolbar.addSpacing(12)
            toolbar.addWidget(QLabel("Axis:", self))
            toolbar.addWidget(self.axis_mode_combo)
            self.reference_lap_combo = QComboBox(self)
            self.reference_lap_combo.setToolTip("Reference lap")
            self.comparison_lap_combo = QComboBox(self)
            self.comparison_lap_combo.setToolTip("Comparison lap")
            self.distance_offset_spin = QDoubleSpinBox(self)
            self.distance_offset_spin.setRange(-10000.0, 10000.0)
            self.distance_offset_spin.setDecimals(3)
            self.distance_offset_spin.setSingleStep(0.5)
            self.distance_offset_spin.setSuffix(" m")
            self.distance_offset_spin.setToolTip("Manual comparison distance offset")
            self.time_offset_spin = QDoubleSpinBox(self)
            self.time_offset_spin.setRange(-10000.0, 10000.0)
            self.time_offset_spin.setDecimals(3)
            self.time_offset_spin.setSingleStep(0.05)
            self.time_offset_spin.setSuffix(" s")
            self.time_offset_spin.setToolTip("Manual comparison time offset")
            toolbar.addWidget(QLabel("Ref:", self))
            toolbar.addWidget(self.reference_lap_combo)
            toolbar.addWidget(QLabel("Compare:", self))
            toolbar.addWidget(self.comparison_lap_combo)
            toolbar.addWidget(QLabel("Δs:", self))
            toolbar.addWidget(self.distance_offset_spin)
            toolbar.addWidget(QLabel("Δt:", self))
            toolbar.addWidget(self.time_offset_spin)
            toolbar.addWidget(self.play_button)
            toolbar.addWidget(self.speed_combo)
            toolbar.addStretch(1)
            toolbar.addWidget(self.cursor_time_label)
            toolbar.addSpacing(12)
            toolbar.addWidget(self.delta_label)
            root.addLayout(toolbar)

            self.case_lap_tree = QTreeWidget(self)
            self.case_lap_tree.setHeaderLabels(["Cases and laps"])
            self.case_lap_tree.setMinimumWidth(220)
            self.channel_selector = ChannelSelector(self)
            browser_layout = QVBoxLayout()
            browser_layout.setContentsMargins(0, 0, 0, 0)
            browser_layout.addWidget(QLabel("Result browser", self))
            browser_layout.addWidget(self.case_lap_tree, 1)
            browser_layout.addWidget(QLabel("Channels", self))
            browser_layout.addWidget(self.channel_selector, 2)
            browser = QWidget(self)
            browser.setLayout(browser_layout)

            self.plot_stack = PlotStack(self)
            self.map_widget = TrackMapWidget(self)
            self.quality_panel = QualityPanel(self)
            self.setup_table = SetupComparisonTable(self)
            self.status_label = QLabel("", self)
            self.status_label.setWordWrap(True)
            self.status_label.setStyleSheet("color: #f4c95d;")
            info_tabs = QTabWidget(self)
            info_tabs.addTab(self.quality_panel, "Quality")
            info_tabs.addTab(self.setup_table, "Setup comparison")
            info = QWidget(self)
            info_layout = QVBoxLayout(info)
            info_layout.setContentsMargins(0, 0, 0, 0)
            info_layout.addWidget(self.map_widget, 3)
            info_layout.addWidget(info_tabs, 2)
            info_layout.addWidget(self.status_label)

            splitter = QSplitter(Qt.Orientation.Horizontal, self)
            splitter.addWidget(browser)
            splitter.addWidget(self.plot_stack)
            splitter.addWidget(info)
            splitter.setStretchFactor(0, 0)
            splitter.setStretchFactor(1, 1)
            splitter.setStretchFactor(2, 0)
            splitter.setSizes([260, 800, 380])
            root.addWidget(splitter, 1)

            self.cursor_slider = QSlider(Qt.Orientation.Horizontal, self)
            self.cursor_slider.setObjectName("cursorSlider")
            self.cursor_slider.setMinimum(0)
            self.cursor_slider.setMaximum(0)
            root.addWidget(self.cursor_slider)

            self.playback_controller = PlaybackController(0, self)
            self.open_button.clicked.connect(self.open_result_dialog)
            self.export_csv_button.clicked.connect(self.export_csv_dialog)
            self.export_screenshot_button.clicked.connect(self.export_screenshot_dialog)
            self.axis_mode_combo.currentIndexChanged.connect(self._axis_changed)
            self.reference_lap_combo.currentIndexChanged.connect(self._reference_lap_changed)
            self.comparison_lap_combo.currentIndexChanged.connect(self._comparison_lap_changed)
            self.distance_offset_spin.valueChanged.connect(self._offsets_changed)
            self.time_offset_spin.valueChanged.connect(self._offsets_changed)
            self.play_button.clicked.connect(self.playback_controller.toggle)
            self.speed_combo.currentIndexChanged.connect(self._speed_changed)
            self.cursor_slider.valueChanged.connect(self.set_cursor_index)
            self.case_lap_tree.itemClicked.connect(self._browser_item_clicked)
            self.channel_selector.channelsChanged.connect(self._channels_changed)
            self.playback_controller.indexChanged.connect(self._playback_index_changed)
            self.playback_controller.playingChanged.connect(self._playing_changed)

        def _populate_browser(self) -> None:
            self.case_lap_tree.clear()
            multiple_results = len(self.model.results) > 1
            for result_name, result in self.model.results.items():
                result_parent: QTreeWidgetItem | QTreeWidget = self.case_lap_tree
                if multiple_results:
                    result_parent = QTreeWidgetItem(self.case_lap_tree, [result_name])
                    result_parent.setData(0, Qt.ItemDataRole.UserRole, (result_name, None, None))
                for case_id, case in result.cases.items():
                    case_label = str(case.metadata.get("label", case_id))
                    case_item = QTreeWidgetItem(result_parent, [case_label])
                    case_item.setData(0, Qt.ItemDataRole.UserRole, (result_name, case_id, None))
                    case_item.setExpanded(True)
                    for lap_id, lap in case.laps.items():
                        role = lap.metadata.get("role", lap_id)
                        lap_item = QTreeWidgetItem(case_item, [f"{role} · {lap_id}"])
                        lap_item.setData(0, Qt.ItemDataRole.UserRole, (result_name, case_id, lap_id))
            self.case_lap_tree.expandAll()

        def _populate_lap_combos(self) -> None:
            case = self.model.current_case
            self.reference_lap_combo.blockSignals(True)
            self.comparison_lap_combo.blockSignals(True)
            try:
                self.reference_lap_combo.clear()
                self.comparison_lap_combo.clear()
                self.comparison_lap_combo.addItem("(none)", None)
                if case is None:
                    return
                for lap_id, lap in case.laps.items():
                    label = f"{lap.metadata.get('role', lap_id)} · {lap_id}"
                    self.reference_lap_combo.addItem(label, lap_id)
                    self.comparison_lap_combo.addItem(label, lap_id)
                reference_index = self.reference_lap_combo.findData(self.model.reference_lap_id)
                comparison_index = self.comparison_lap_combo.findData(self.model.comparison_lap_id)
                self.reference_lap_combo.setCurrentIndex(max(0, reference_index))
                self.comparison_lap_combo.setCurrentIndex(max(0, comparison_index))
            finally:
                self.reference_lap_combo.blockSignals(False)
                self.comparison_lap_combo.blockSignals(False)

        def _refresh_view(self, *, reset_playback: bool = False) -> None:
            self._populate_lap_combos()
            lap = self.model.reference_lap or self.model.current_lap
            if reset_playback:
                self.playback_controller.pause()
            sample_count = len(_preferred_axis(lap).time_s) if _preferred_axis(lap) is not None else 0
            self.playback_controller.set_sample_count(sample_count, reset=reset_playback)
            self.cursor_slider.setMaximum(max(0, sample_count - 1))
            self.channel_selector.set_channels(self.model.available_channels, self.model.selected_channel_ids)
            self.map_widget.set_track(self.model.track)
            self.plot_stack.set_data(
                self.model.reference_lap or self.model.current_lap,
                self.model.selected_channel_ids,
                comparison_lap=self.model.comparison_lap,
                axis_mode=self.model.axis_mode,
                cursor_value=self.model.cursor_value,
                comparison_offset=self.model.distance_offset_m if self.model.axis_mode == "distance" else self.model.time_offset_s,
            )
            self._refresh_info()
            self._update_cursor_widgets()
            if lap is None:
                self.status_label.setText("Open a QSS telemetry result or inject a ResultFile to begin.")
            else:
                self._update_status()

        def _refresh_info(self) -> None:
            self.quality_panel.update_quality(
                self.model.reference_lap or self.model.current_lap,
                self.model.track,
                missing_channels=self.model.missing_channels,
                comparison_status=self.model.comparison_status,
                distance_offset_m=self.model.distance_offset_m,
                time_offset_s=self.model.time_offset_s,
            )
            reference_setup, comparison_setup = self.model.setup_comparison
            self.setup_table.set_setups(reference_setup, comparison_setup)
            self.delta_label.set_delta(self.model.delta_seconds)

        def _update_status(self) -> None:
            if self.model.missing_channels:
                self.status_label.setText("Selected channel(s) not available: " + ", ".join(self.model.missing_channels))
            elif self.model.comparison_lap is not None and self.model.comparison_error:
                self.status_label.setText("Comparison unavailable: " + self.model.comparison_error)
            elif self.model.comparison_lap is not None:
                self.status_label.setText("")

        def _update_cursor_widgets(self) -> None:
            self.cursor_slider.blockSignals(True)
            self.cursor_slider.setValue(self.model.cursor_index)
            self.cursor_slider.blockSignals(False)
            time_s = self.model.cursor_time_s
            distance_m = self.model.cursor_distance_m
            time_text = "—" if time_s is None else f"{time_s:.3f} s"
            distance_text = "—" if distance_m is None else f"{distance_m:.3f} m"
            self.cursor_time_label.setText(f"t={time_text} · s={distance_text}")
            self.plot_stack.set_cursor_value(self.model.cursor_value)
            self.map_widget.set_marker_distance(distance_m)

        def _browser_item_clicked(self, item: QTreeWidgetItem, column: int) -> None:
            del column
            data = item.data(0, Qt.ItemDataRole.UserRole)
            if not data or data[1] is None or data[2] is None:
                return
            result_name, case_id, lap_id = data
            if result_name != self.model.result_name:
                self.model.select_result(result_name)
            self.model.select_lap(case_id, lap_id)
            self._refresh_view(reset_playback=True)

        def _channels_changed(self, channel_ids: list[str]) -> None:
            self.model.set_selected_channels(channel_ids)
            self.plot_stack.set_data(
                self.model.reference_lap or self.model.current_lap,
                self.model.selected_channel_ids,
                comparison_lap=self.model.comparison_lap,
                axis_mode=self.model.axis_mode,
                cursor_value=self.model.cursor_value,
                comparison_offset=self.model.distance_offset_m if self.model.axis_mode == "distance" else self.model.time_offset_s,
            )
            self._refresh_info()
            self._update_status()

        def _playback_index_changed(self, index: int) -> None:
            self.model.set_cursor_index(index)
            self._update_cursor_widgets()

        def _playing_changed(self, playing: bool) -> None:
            self.play_button.setText("Pause" if playing else "Play")

        def _speed_changed(self, index: int) -> None:
            value = self.speed_combo.itemData(index)
            if value is not None:
                self.playback_controller.set_speed(float(value))

        def _offsets_changed(self, value: float) -> None:
            del value
            self.model.set_alignment_offsets(
                distance_offset_m=self.distance_offset_spin.value(),
                time_offset_s=self.time_offset_spin.value(),
            )
            self.plot_stack.set_data(
                self.model.reference_lap or self.model.current_lap,
                self.model.selected_channel_ids,
                comparison_lap=self.model.comparison_lap,
                axis_mode=self.model.axis_mode,
                cursor_value=self.model.cursor_value,
                comparison_offset=self.model.distance_offset_m if self.model.axis_mode == "distance" else self.model.time_offset_s,
            )
            self._refresh_info()
            self._update_status()

        def _reference_lap_changed(self, index: int) -> None:
            if index < 0:
                return
            reference_id = self.reference_lap_combo.itemData(index)
            if reference_id is None:
                return
            comparison_id = self.comparison_lap_combo.currentData()
            self.model.set_comparison_laps(str(reference_id), str(comparison_id) if comparison_id is not None else None)
            self._refresh_view(reset_playback=True)

        def _comparison_lap_changed(self, index: int) -> None:
            reference_id = self.reference_lap_combo.currentData()
            if reference_id is None:
                return
            comparison_id = self.comparison_lap_combo.itemData(index)
            self.model.set_comparison_laps(str(reference_id), str(comparison_id) if comparison_id is not None else None)
            self._refresh_view(reset_playback=False)

        def _axis_changed(self, index: int) -> None:
            mode = self.axis_mode_combo.itemData(index)
            if mode is None:
                return
            self.model.set_axis_mode(str(mode))
            self.plot_stack.set_data(
                self.model.reference_lap or self.model.current_lap,
                self.model.selected_channel_ids,
                comparison_lap=self.model.comparison_lap,
                axis_mode=self.model.axis_mode,
                cursor_value=self.model.cursor_value,
                comparison_offset=self.model.distance_offset_m if self.model.axis_mode == "distance" else self.model.time_offset_s,
            )
            self._update_cursor_widgets()
            self._update_status()

        def set_cursor_index(self, index: int) -> int:
            bounded = self.model.set_cursor_index(index)
            self.playback_controller.set_index(bounded)
            self._update_cursor_widgets()
            return bounded

        def set_cursor_value(self, value: float) -> int:
            index = self.model.set_cursor_value(value)
            self.playback_controller.set_index(index)
            self._update_cursor_widgets()
            return index

        def set_axis_mode(self, mode: str) -> None:
            index = self.axis_mode_combo.findData("distance" if mode == "distance" else "time")
            if index >= 0:
                self.axis_mode_combo.setCurrentIndex(index)
            else:
                self.model.set_axis_mode(mode)
                self._axis_changed(-1)

        def set_selected_channels(self, channel_ids: Iterable[str]) -> list[str]:
            requested = list(channel_ids)
            self.model.set_selected_channels(requested)
            self.channel_selector.set_selected(self.model.selected_channel_ids, emit=False)
            self._channels_changed(requested)
            return self.model.missing_channels

        def set_comparison_laps(self, reference_lap_id: str, comparison_lap_id: str | None, case_id: str | None = None) -> None:
            self.model.set_comparison_laps(reference_lap_id, comparison_lap_id, case_id)
            self._refresh_view(reset_playback=False)

        def set_alignment_offsets(self, *, distance_offset_m: float | None = None, time_offset_s: float | None = None) -> None:
            if distance_offset_m is not None:
                self.distance_offset_spin.setValue(float(distance_offset_m))
            if time_offset_s is not None:
                self.time_offset_spin.setValue(float(time_offset_s))
            self.model.set_alignment_offsets(
                distance_offset_m=self.distance_offset_spin.value(),
                time_offset_s=self.time_offset_spin.value(),
            )
            self._offsets_changed(0.0)

        def open_result_dialog(self) -> None:
            path, _ = QFileDialog.getOpenFileName(self, "Open QSS telemetry result", "", "QSS HDF5 (*.h5 *.hdf5 *.qss.h5);;All files (*)")
            if path:
                self.load_result_file(path)

        def load_result_file(self, path: str | Path) -> None:
            try:
                from .io import load_results  # lazy core adapter

                loaded = load_results(path)
                incoming = _result_map(loaded)
                self.model.results.update(incoming)
                if incoming:
                    self.model.result_name = next(reversed(incoming))
                    self.model.case_id, self.model.lap_id = _first_case_lap(self.model.current_result)
                    self.model.reference_lap_id = self.model.lap_id
                    self.model.comparison_lap_id = None
                    self.model._select_default_channels()
                self._populate_browser()
                self._refresh_view(reset_playback=True)
                self.status_label.setText(f"Loaded {Path(path).name}")
            except Exception as error:  # surface reader/validation errors in the local UI
                self.status_label.setText(f"Could not load {Path(path).name}: {error}")

        def export_csv_dialog(self) -> None:
            path, _ = QFileDialog.getSaveFileName(self, "Export selected telemetry", "telemetry.csv", "CSV (*.csv)")
            if path:
                try:
                    self.export_csv(path)
                    self.status_label.setText(f"Exported {Path(path).name}")
                except Exception as error:
                    self.status_label.setText(f"CSV export failed: {error}")

        def export_screenshot_dialog(self) -> None:
            path, _ = QFileDialog.getSaveFileName(self, "Export viewer screenshot", "telemetry.png", "PNG image (*.png);;JPEG image (*.jpg *.jpeg)")
            if path:
                try:
                    self.export_screenshot(path)
                    self.status_label.setText(f"Exported {Path(path).name}")
                except Exception as error:
                    self.status_label.setText(f"Screenshot export failed: {error}")

        def export_csv(self, path: str | Path) -> Path:
            """Export every source sample for the selected channels."""

            lap = self.model.reference_lap or self.model.current_lap
            if lap is None:
                raise ValueError("no lap is selected")
            target = Path(path)
            target.parent.mkdir(parents=True, exist_ok=True)
            axis = _preferred_axis(lap)
            if axis is None:
                raise ValueError("selected lap has no axis")
            selected = [lap.channels[channel_id] for channel_id in self.model.selected_channel_ids if channel_id in lap.channels]
            common = all(channel.values.size == axis.time_s.size and channel.metadata.axis_id == axis.id for channel in selected)
            with target.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                if common:
                    headers = ["time_s", "distance_m"]
                    for channel in selected:
                        headers.extend((channel.metadata.id, f"{channel.metadata.id}__valid"))
                    writer.writerow(headers)
                    for row_index in range(axis.time_s.size):
                        row: list[Any] = [axis.time_s[row_index], axis.distance_m[row_index]]
                        for channel in selected:
                            value = channel.values[row_index] if channel.valid[row_index] else math.nan
                            row.extend((value, int(channel.valid[row_index])))
                        writer.writerow(row)
                else:
                    # Mixed-axis results retain each source channel's full
                    # resolution instead of silently downsampling to one axis.
                    columns = ["time_s", "distance_m"]
                    for channel in selected:
                        columns.extend((f"{channel.metadata.id}__time_s", f"{channel.metadata.id}__distance_m", channel.metadata.id, f"{channel.metadata.id}__valid"))
                    writer.writerow(columns)
                    row_count = max([axis.time_s.size] + [channel.values.size for channel in selected])
                    for row_index in range(row_count):
                        row = [axis.time_s[row_index] if row_index < axis.time_s.size else math.nan, axis.distance_m[row_index] if row_index < axis.distance_m.size else math.nan]
                        for channel in selected:
                            channel_axis = lap.axes[channel.metadata.axis_id]
                            if row_index < channel.values.size:
                                row.extend((channel_axis.time_s[row_index], channel_axis.distance_m[row_index], channel.values[row_index] if channel.valid[row_index] else math.nan, int(channel.valid[row_index])))
                            else:
                                row.extend((math.nan, math.nan, math.nan, 0))
                        writer.writerow(row)
            return target

        def export_screenshot(self, path: str | Path) -> Path:
            target = Path(path)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not self.grab().save(str(target)):
                raise OSError(f"could not save screenshot to {target}")
            return target

        def closeEvent(self, event: Any) -> None:
            self.playback_controller.pause()
            super().closeEvent(event)


    # Friendly aliases for callers that prefer a shorter application name.
    TelemetryViewer = TelemetryViewerWindow
    ViewerWindow = TelemetryViewerWindow

else:

    class TelemetryViewerWindow:  # pragma: no cover - dependency guard
        """Dependency guard that leaves the data model usable without Qt."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            del args, kwargs
            detail = f": {_QT_IMPORT_ERROR}" if _QT_IMPORT_ERROR is not None else ""
            raise RuntimeError("PySide6 and pyqtgraph are required for the desktop viewer" + detail)


    TelemetryViewer = TelemetryViewerWindow
    ViewerWindow = TelemetryViewerWindow


def main(argv: Sequence[str] | None = None) -> int:
    """Launch the viewer and lazily load a path supplied on the command line."""

    if not _QT_AVAILABLE:
        detail = f": {_QT_IMPORT_ERROR}" if _QT_IMPORT_ERROR is not None else ""
        raise RuntimeError("PySide6 and pyqtgraph are required for qss-telemetry-viewer" + detail)
    arguments = list(sys.argv[1:] if argv is None else argv)
    application = QApplication.instance() or QApplication([sys.argv[0], *arguments])
    window = TelemetryViewerWindow()
    window.show()
    if arguments:
        window.load_result_file(arguments[0])
    return application.exec()


__all__ = [
    "TelemetryViewer",
    "TelemetryViewerWindow",
    "ViewerModel",
    "ViewerWindow",
    "delta_time_seconds",
    "integrate_curvature",
    "main",
    "track_geometry",
]
