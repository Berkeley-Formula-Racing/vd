"""Grouped, searchable channel selection widget."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QLineEdit, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget

from ..schema import ChannelData


def channel_group(channel_id: str, channel: ChannelData | None = None) -> str:
    """Return a friendly group for a self-describing channel."""

    metadata = channel.metadata if channel is not None else None
    explicit = ""
    if metadata is not None:
        # Exporters may provide a group as an additional metadata field.  The
        # shared v1 fields remain sufficient when it is absent.
        explicit = str(getattr(metadata, "group", "") or "")
        if not explicit and hasattr(metadata, "as_dict"):
            explicit = str(metadata.as_dict().get("group", "") or "")
    if explicit:
        return explicit.title()

    value = channel_id.lower()
    label = metadata.label.lower() if metadata is not None else ""
    haystack = f"{value} {label}"
    if any(token in haystack for token in ("speed", "velocity", "accel", "yaw", "slip", "kinematic")):
        return "Vehicle dynamics"
    if any(token in haystack for token in ("steer", "throttle", "brake", "control", "demand", "pedal", "torque")):
        return "Controls"
    if any(token in haystack for token in ("wheel", "tire", "tyre", "suspension", "damper", "load")):
        return "Wheels and suspension"
    if any(token in haystack for token in ("delta", "lap", "sector", "time")):
        return "Lap metrics"
    return "Other"


class ChannelSelector(QWidget):
    """A checkable channel tree with a search field and group checkboxes.

    The widget emits only channel IDs, so the rest of the UI stays independent
    of the Qt item hierarchy and can be driven directly in tests.
    """

    channelsChanged = Signal(list)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.search_edit = QLineEdit(self)
        self.search_edit.setPlaceholderText("Search channels")
        self.search_edit.setClearButtonEnabled(True)
        self.search_edit.textChanged.connect(self._filter_items)

        self.tree = QTreeWidget(self)
        self.tree.setHeaderLabels(["Channel", "Unit"])
        self.tree.setRootIsDecorated(True)
        self.tree.setAlternatingRowColors(True)
        self.tree.itemChanged.connect(self._item_changed)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.search_edit)
        layout.addWidget(self.tree, 1)

        self._channels: dict[str, ChannelData] = {}
        self._channel_items: dict[str, QTreeWidgetItem] = {}
        self._group_items: dict[str, QTreeWidgetItem] = {}
        self._updating = False

    def set_channels(
        self,
        channels: Mapping[str, ChannelData] | Iterable[ChannelData],
        selected: Iterable[str] | None = None,
    ) -> None:
        if isinstance(channels, Mapping):
            values = list(channels.values())
        else:
            values = list(channels)
        self._channels = {item.metadata.id: item for item in values}
        selected_set = set(self._channels) if selected is None else set(selected) & set(self._channels)
        self._updating = True
        try:
            self.tree.clear()
            self._channel_items.clear()
            self._group_items.clear()
            groups: dict[str, list[ChannelData]] = {}
            for channel in sorted(values, key=lambda item: (channel_group(item.metadata.id, item), item.metadata.label.lower(), item.metadata.id)):
                groups.setdefault(channel_group(channel.metadata.id, channel), []).append(channel)
            for group_name, group_channels in groups.items():
                group_item = QTreeWidgetItem(self.tree, [group_name, ""])
                group_item.setFlags(group_item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                group_item.setCheckState(0, Qt.CheckState.Unchecked)
                self._group_items[group_name] = group_item
                for channel in group_channels:
                    metadata = channel.metadata
                    item = QTreeWidgetItem(group_item, [metadata.label or metadata.id, metadata.unit])
                    item.setData(0, Qt.ItemDataRole.UserRole, metadata.id)
                    item.setToolTip(0, metadata.description or metadata.id)
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(0, Qt.CheckState.Checked if metadata.id in selected_set else Qt.CheckState.Unchecked)
                    self._channel_items[metadata.id] = item
                self._update_group_state(group_item)
        finally:
            self._updating = False
        self._filter_items(self.search_edit.text())

    def checked_channels(self) -> list[str]:
        return [
            channel_id
            for channel_id, item in self._channel_items.items()
            if item.checkState(0) == Qt.CheckState.Checked
        ]

    def is_checked(self, channel_id: str) -> bool:
        item = self._channel_items.get(channel_id)
        return item is not None and item.checkState(0) == Qt.CheckState.Checked

    def set_checked(self, channel_id: str, checked: bool, *, emit: bool = True) -> None:
        item = self._channel_items.get(channel_id)
        if item is None:
            return
        self._updating = True
        try:
            item.setCheckState(0, Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked)
            self._update_group_state(item.parent())
        finally:
            self._updating = False
        if emit:
            self.channelsChanged.emit(self.checked_channels())

    def set_selected(self, channel_ids: Iterable[str], *, emit: bool = True) -> None:
        selected = set(channel_ids)
        self._updating = True
        try:
            for channel_id, item in self._channel_items.items():
                item.setCheckState(0, Qt.CheckState.Checked if channel_id in selected else Qt.CheckState.Unchecked)
            for group in self._group_items.values():
                self._update_group_state(group)
        finally:
            self._updating = False
        if emit:
            self.channelsChanged.emit(self.checked_channels())

    def channel_ids(self) -> list[str]:
        return list(self._channels)

    def _item_changed(self, item: QTreeWidgetItem, column: int) -> None:
        if self._updating or column != 0:
            return
        self._updating = True
        try:
            if item in self._group_items.values():
                state = item.checkState(0)
                for index in range(item.childCount()):
                    item.child(index).setCheckState(0, state if state != Qt.CheckState.PartiallyChecked else Qt.CheckState.Unchecked)
            else:
                self._update_group_state(item.parent())
        finally:
            self._updating = False
        self.channelsChanged.emit(self.checked_channels())

    @staticmethod
    def _update_group_state(group: QTreeWidgetItem | None) -> None:
        if group is None or group.childCount() == 0:
            return
        checked = sum(group.child(index).checkState(0) == Qt.CheckState.Checked for index in range(group.childCount()))
        if checked == 0:
            state = Qt.CheckState.Unchecked
        elif checked == group.childCount():
            state = Qt.CheckState.Checked
        else:
            state = Qt.CheckState.PartiallyChecked
        group.setCheckState(0, state)

    def _filter_items(self, text: str) -> None:
        query = text.strip().lower()
        for group_name, group in self._group_items.items():
            visible_children = 0
            for index in range(group.childCount()):
                item = group.child(index)
                channel_id = str(item.data(0, Qt.ItemDataRole.UserRole) or "")
                channel = self._channels.get(channel_id)
                haystack = f"{channel_id} {item.text(0)} {item.text(1)} {(channel.metadata.description if channel else '')}".lower()
                visible = not query or query in haystack
                item.setHidden(not visible)
                visible_children += visible
            group.setHidden(visible_children == 0)
            if visible_children:
                group.setExpanded(True)
