"""Small information panes used beside the telemetry plots."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QAbstractItemView, QLabel, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget

from ..schema import CaseData, LapData, TrackData


def _value_text(value: Any) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (dict, list, tuple)):
        return ", ".join(f"{key}={item}" for key, item in value.items()) if isinstance(value, dict) else str(value)
    return str(value)


class SetupComparisonTable(QTableWidget):
    """Display setup keys side-by-side for the selected lap pair."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setColumnCount(4)
        self.setHorizontalHeaderLabels(["Parameter", "Reference", "Comparison", "Change"])
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setAlternatingRowColors(True)
        self.horizontalHeader().setStretchLastSection(True)

    def set_setups(self, reference: Mapping[str, Any] | None, comparison: Mapping[str, Any] | None) -> None:
        left = dict(reference or {})
        right = dict(comparison or {})
        keys = sorted(set(left) | set(right), key=str.casefold)
        self.setRowCount(len(keys))
        for row, key in enumerate(keys):
            left_text = _value_text(left.get(key))
            right_text = _value_text(right.get(key))
            changed = "" if left_text == right_text else "changed"
            for column, text in enumerate((str(key), left_text, right_text, changed)):
                self.setItem(row, column, QTableWidgetItem(text))
            if changed:
                self.item(row, 3).setForeground(Qt.GlobalColor.yellow)
        self.resizeColumnsToContents()


class QualityPanel(QWidget):
    """Show validity, reconstruction, missing-channel, and map quality facts."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.map_label = QLabel("Map: —", self)
        self.comparison_label = QLabel("Comparison: —", self)
        self.reconstruction_label = QLabel("Reconstruction: —", self)
        self.missing_label = QLabel("Missing channels: none", self)
        self.validity_label = QLabel("Validity: —", self)
        self.notes_label = QLabel("", self)
        self.notes_label.setWordWrap(True)
        self.notes_label.setStyleSheet("color: #9aa9bd;")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        for widget in (self.map_label, self.comparison_label, self.reconstruction_label, self.missing_label, self.validity_label, self.notes_label):
            layout.addWidget(widget)
        layout.addStretch(1)

    def update_quality(
        self,
        lap: LapData | None,
        track: TrackData | None,
        *,
        missing_channels: Iterable[str] = (),
        comparison_status: str | None = None,
        distance_offset_m: float = 0.0,
        time_offset_s: float = 0.0,
    ) -> None:
        missing = list(missing_channels)
        if track is None:
            map_quality = "—"
        elif track.x_m is not None and track.y_m is not None:
            map_quality = "Supplied geometry"
        else:
            map_quality = "Schematic map"
        self.map_label.setText(f"Map: {map_quality}")
        status = comparison_status or "—"
        self.comparison_label.setText(f"Comparison: {status} (Δs={distance_offset_m:.3f} m, Δt={time_offset_s:.3f} s)")
        metadata = dict(lap.metadata) if lap is not None else {}
        reconstruction = metadata.get("reconstruction") or metadata.get("qss_reconstruction") or metadata.get("reconstruction_status")
        self.reconstruction_label.setText(f"Reconstruction: {_value_text(reconstruction) if reconstruction is not None else 'not reported'}")
        self.missing_label.setText("Missing channels: " + (", ".join(missing) if missing else "none"))
        if lap is None:
            self.validity_label.setText("Validity: —")
            self.notes_label.setText("")
            return
        total = 0
        invalid = 0
        for channel in lap.channels.values():
            total += channel.valid.size
            invalid += int((~channel.valid).sum())
        valid_percent = 100.0 if total == 0 else 100.0 * (total - invalid) / total
        self.validity_label.setText(f"Validity: {valid_percent:.1f}% ({invalid} invalid samples)")
        notes = metadata.get("failure_reason") or metadata.get("quality_note") or ""
        self.notes_label.setText(str(notes))


class DeltaDisplay(QLabel):
    """A label with a shared sign convention: positive means slower."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__("Delta: —", parent)
        self.setObjectName("deltaDisplay")

    def set_delta(self, seconds: float | None) -> None:
        if seconds is None:
            self.setText("Delta: —")
            return
        value = float(seconds)
        if abs(value) < 5e-4:
            self.setText("Delta: 0.000 s (equal)")
        elif value > 0:
            self.setText(f"Delta: +{value:.3f} s (comparison slower)")
        else:
            self.setText(f"Delta: {value:.3f} s (comparison faster)")
