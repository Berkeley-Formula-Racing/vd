"""Tkinter tree and record editor for the offline vehicle test log."""

from __future__ import annotations

import os
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from .event_ui import EventLayoutManager
from .models import EventLayout, Lap, Setup, StagedAttachment, TestDay, new_id, utc_now_iso
from .time_value import format_lap_time, parse_lap_time
from .validation import ValidationError


class TestLogWindow:
    """Main three-level test-day, setup, and lap editor."""

    def __init__(self, root: tk.Tk, services):
        self.root = root
        self.services = services
        self.tree: ttk.Treeview | None = None
        self.fields: dict[str, tk.Widget] = {}
        self.actions: dict[str, ttk.Button] = {}
        self.attachment_list: tk.Listbox | None = None
        self.status_label: ttk.Label | None = None
        self.dirty_label: ttk.Label | None = None
        self.event_manager: EventLayoutManager | None = None
        self._editor_frame: ttk.Frame | None = None
        self._current_kind: str | None = None
        self._record: TestDay | Setup | Lap | None = None
        self._is_new = False
        self._context_id: str | None = None
        self._loaded_values: dict[str, str] = {}
        self._staged: list[StagedAttachment] = []
        self._attachment_refs: list[tuple[str, object]] = []
        self._field_errors: dict[str, ttk.Label] = {}
        self._event_id_by_label: dict[str, str] = {}
        self._event_label_by_id: dict[str, str] = {}
        self._lap_order_original: list[str] = []
        self._pending_lap_order: list[str] | None = None
        self._suppress_tree_selection = False
        self._built = False
        self._closed = False

    @property
    def staged_attachments(self) -> list[StagedAttachment]:
        """Staged file copies currently held by the editor."""
        return list(self._staged)

    def build(self) -> None:
        if self._built:
            return
        self._built = True
        self.root.title("VD Vehicle Test Log")
        self.root.geometry("1180x720")
        self.root.minsize(900, 630)
        self.root.protocol("WM_DELETE_WINDOW", self.close_request)

        toolbar = ttk.Frame(self.root, padding=(8, 8, 8, 4))
        toolbar.pack(fill="x")
        button_specs = (
            ("add_day", "Add Test Day", self.add_day),
            ("add_setup", "Add Setup", self.add_setup),
            ("add_lap", "Add Lap", self.add_lap),
            ("duplicate_setup", "Duplicate Setup", self.duplicate_setup),
            ("manage_events", "Manage Event/Layouts", self.manage_events),
            ("attach_file", "Attach File", self.attach_file),
            ("open_attachment", "Open Attachment", self.open_attachment),
            ("open_data_folder", "Open Data Folder", self.open_data_folder),
            ("move_lap_up", "Move Lap Up", lambda: self.move_lap(-1)),
            ("move_lap_down", "Move Lap Down", lambda: self.move_lap(1)),
            ("save", "Save", self.save_current),
            ("delete", "Delete", self.delete_current),
            ("backup", "Backup", self.create_backup),
            ("export_csv", "Export CSV", self.export_csv),
        )
        for index, (name, label, callback) in enumerate(button_specs):
            row, column = divmod(index, 7)
            button = ttk.Button(toolbar, text=label, command=callback)
            button.grid(row=row, column=column, sticky="ew", padx=(0, 5), pady=(0, 4))
            self.actions[name] = button
        for column in range(7):
            toolbar.columnconfigure(column, weight=1)

        body = ttk.Panedwindow(self.root, orient="horizontal")
        body.pack(fill="both", expand=True, padx=8, pady=4)
        tree_frame = ttk.Frame(body, padding=(0, 0, 8, 0))
        self._editor_frame = ttk.Frame(body, padding=(8, 0, 0, 0))
        body.add(tree_frame, weight=1)
        body.add(self._editor_frame, weight=2)

        self.tree = ttk.Treeview(tree_frame, show="tree", selectmode="browse")
        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.tree.bind("<<TreeviewSelect>>", self._on_tree_select)

        bottom = ttk.Frame(self.root, padding=(8, 4, 8, 8))
        bottom.pack(fill="x")
        self.dirty_label = ttk.Label(bottom, text="")
        self.dirty_label.pack(side="left", padx=(0, 12))
        self.status_label = ttk.Label(bottom, text="", anchor="w", wraplength=950)
        self.status_label.pack(side="left", fill="x", expand=True)

        self._show_editor(None, None, is_new=True)
        self.refresh_tree()
        self._update_action_states()

    def _set_status(self, text: str) -> None:
        if self.status_label is not None:
            self.status_label.configure(text=text)

    def _set_dirty_text(self, text: str) -> None:
        if self.dirty_label is not None:
            self.dirty_label.configure(text=text)


    def _add_form_field(self, parent, key: str, label: str, row: int, *, multiline=False, choices=None) -> int:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="nw" if multiline else "w", padx=(0, 8), pady=5)
        if choices is not None:
            widget = ttk.Combobox(parent, values=choices, state="readonly")
        elif multiline:
            widget = tk.Text(parent, height=4 if key != "settings_text" else 8, width=38, wrap="word")
        else:
            widget = ttk.Entry(parent)
        widget.grid(row=row, column=1, sticky="ew", pady=5)
        self.fields[key] = widget
        error_label = ttk.Label(parent, text="", foreground="#a02020", wraplength=280)
        error_label.grid(row=row + 1, column=1, sticky="w")
        self._field_errors[key] = error_label
        if isinstance(widget, ttk.Combobox):
            widget.bind("<<ComboboxSelected>>", self._on_field_change)
        elif isinstance(widget, tk.Text):
            widget.bind("<<Modified>>", self._on_text_modified)
        else:
            widget.bind("<KeyRelease>", self._on_field_change)
        return row + 2

    def _show_editor(self, kind: str | None, record, *, is_new: bool, context_id: str | None = None) -> None:
        if self._editor_frame is None:
            return
        for child in self._editor_frame.winfo_children():
            child.destroy()
        self.fields = {}
        self._field_errors = {}
        self._attachment_refs = []
        self._current_kind = kind
        self._record = record
        self._is_new = is_new
        self._context_id = context_id
        self._lap_order_original = []
        self._pending_lap_order = None
        if kind == "lap":
            setup_id = record.setup_id if isinstance(record, Lap) else context_id
            if setup_id:
                self._lap_order_original = [lap.id for lap in self.services.list_laps(setup_id)]
                self._pending_lap_order = list(self._lap_order_original)
        self._staged.clear()
        self._event_id_by_label.clear()
        self._event_label_by_id.clear()
        self._set_status("")
        if kind is None:
            self._loaded_values = {}
            self.attachment_list = None
            self._set_dirty_text("")
            self._update_action_states()
            return

        panel = ttk.Frame(self._editor_frame)
        panel.pack(fill="both", expand=True)
        panel.columnconfigure(0, weight=1)
        panel.rowconfigure(1, weight=1)
        title = {"day": "Test Day", "setup": "Vehicle Setup", "lap": "Lap Record"}[kind]
        ttk.Label(panel, text=("New " if is_new else "Edit ") + title, font=("Segoe UI", 12, "bold")).grid(
            row=0, column=0, sticky="w", pady=(0, 6)
        )
        body = ttk.Frame(panel)
        body.grid(row=1, column=0, sticky="nsew")
        body.columnconfigure(1, weight=1)
        row = 0
        if kind == "day":
            values = self._day_values(record) if record is not None else {}
            row = self._add_form_field(body, "date", "Date (YYYY-MM-DD)", row)
            row = self._add_form_field(body, "location", "Location", row)
            row = self._add_form_field(body, "weather", "Weather", row)
            row = self._add_form_field(body, "notes", "Notes", row, multiline=True)
        elif kind == "setup":
            values = self._setup_values(record) if record is not None else {}
            row = self._add_form_field(body, "name", "Setup label", row)
            row = self._add_form_field(body, "setup_code", "Setup ID", row)
            row = self._add_form_field(body, "settings_text", "Settings", row, multiline=True)
            row = self._add_form_field(body, "notes", "Notes", row, multiline=True)
        else:
            values = self._lap_values(record) if record is not None else {}
            row = self._add_form_field(body, "lap_time", "Lap time (seconds or m:ss.sss)", row)
            row = self._add_form_field(body, "status", "Status", row, choices=("valid", "invalid"))
            row = self._add_form_field(body, "driver", "Driver", row)
            row = self._add_form_field(body, "time_of_day", "Time of day", row)
            row = self._add_form_field(body, "event_layout_id", "Event / layout", row, choices=())
            row = self._add_form_field(body, "notes", "Notes", row, multiline=True)
            self._populate_event_choices(record)
        self._set_form_values(values or {})
        if kind == "lap":
            self._populate_event_choices(record, preserve_value=False)

        attachment_frame = ttk.LabelFrame(panel, text="Attachments", padding=6)
        attachment_frame.grid(row=2, column=0, sticky="nsew", pady=(8, 4))
        self.attachment_list = tk.Listbox(attachment_frame, height=5, exportselection=False)
        self.attachment_list.pack(fill="both", expand=True)
        self._refresh_attachment_list()
        self._loaded_values = self._capture_values()
        self._update_dirty_indicator()
        self._update_action_states()

    def _day_values(self, record: TestDay) -> dict[str, str]:
        return {
            "date": record.date,
            "location": record.location,
            "weather": record.weather or "",
            "notes": record.notes or "",
        }

    def _setup_values(self, record: Setup) -> dict[str, str]:
        return {
            "name": record.name,
            "setup_code": record.setup_code or "",
            "settings_text": record.settings_text,
            "notes": record.notes or "",
        }

    def _lap_values(self, record: Lap) -> dict[str, str]:
        return {
            "lap_time": format_lap_time(record.time_ms),
            "status": record.status,
            "driver": record.driver or "",
            "time_of_day": record.time_of_day or "",
            "event_layout_id": "",
            "notes": record.notes or "",
        }

    def _set_form_values(self, values: dict[str, str]) -> None:
        for key, widget in self.fields.items():
            value = values.get(key, "")
            if isinstance(widget, tk.Text):
                widget.delete("1.0", "end")
                widget.insert("1.0", value)
                widget.edit_modified(False)
            elif isinstance(widget, ttk.Combobox):
                widget.set(value)
            else:
                widget.delete(0, "end")
                widget.insert(0, value)
        self._clear_field_errors()

    def _capture_values(self) -> dict[str, str]:
        values: dict[str, str] = {}
        for key, widget in self.fields.items():
            values[key] = widget.get("1.0", "end-1c") if isinstance(widget, tk.Text) else widget.get()
        return values

    def _clear_field_errors(self) -> None:
        for label in self._field_errors.values():
            label.configure(text="")

    def _show_validation_error(self, error: ValidationError) -> None:
        field = "lap_time" if error.field in ("time_ms", "lap_time") else error.field
        label = self._field_errors.get(field)
        if label is not None:
            label.configure(text=error.message)
        self._set_status(error.message)


    def _on_field_change(self, _event=None) -> None:
        self._update_dirty_indicator()
        self._update_action_states()

    def _on_text_modified(self, event) -> None:
        widget = event.widget
        if widget.edit_modified():
            widget.edit_modified(False)
            self._update_dirty_indicator()
            self._update_action_states()

    def _optional(self, value: str) -> str | None:
        value = value.strip()
        return value or None

    def _populate_event_choices(self, record: Lap | None, *, preserve_value=True) -> None:
        widget = self.fields.get("event_layout_id")
        if not isinstance(widget, ttk.Combobox):
            return
        current_value = widget.get() if preserve_value else ""
        events = self.services.list_event_layouts(include_archived=False)
        if record is not None:
            selected = self.services.get_event_layout(record.event_layout_id)
            if selected is not None and selected.archived and selected not in events:
                events.append(selected)
        displays = [
            (f"{event.event_name + ' · ' if event.event_name else ''}{event.track_name} / {event.layout_name}", event.id)
            for event in events
        ]
        counts: dict[str, int] = {}
        for display, _identifier in displays:
            counts[display] = counts.get(display, 0) + 1
        self._event_id_by_label.clear()
        self._event_label_by_id.clear()
        labels = []
        for display, identifier in displays:
            label = f"{display} [{identifier}]" if counts[display] > 1 else display
            labels.append(label)
            self._event_id_by_label[label] = identifier
            self._event_label_by_id[identifier] = label
        widget.configure(values=labels)
        if current_value:
            widget.set(current_value)
        elif record is not None:
            widget.set(self._event_label_by_id.get(record.event_layout_id, ""))
        else:
            widget.set("")

    def event_choice_label(self, identifier: str) -> str | None:
        """Return the current human-readable choice for an event/layout ID."""
        return self._event_label_by_id.get(identifier)

    def _refresh_event_choices(self) -> None:
        if self._current_kind == "lap":
            record = self._record if isinstance(self._record, Lap) else None
            self._populate_event_choices(record)

    def _tree_item_for_current(self) -> str | None:
        if self._current_kind is None:
            return None
        if not self._is_new and self._record is not None:
            return f"{self._current_kind}:{self._record.id}"
        if self._current_kind == "setup":
            return f"day:{self._context_id}" if self._context_id else None
        if self._current_kind == "lap":
            return f"setup:{self._context_id}" if self._context_id else None
        return None

    def refresh_tree(self, select_item: str | None = None) -> None:
        if self.tree is None:
            return
        old_selection = self.tree.selection()
        desired = select_item or (old_selection[0] if old_selection else self._tree_item_for_current())
        self._suppress_tree_selection = True
        for item in self.tree.get_children(""):
            self.tree.delete(item)
        days = sorted(self.services.list_days(), key=lambda day: (day.date, day.created_at), reverse=True)
        for day in days:
            day_item = f"day:{day.id}"
            self.tree.insert("", "end", iid=day_item, text=f"{day.date} · {day.location}", open=True)
            for setup in sorted(self.services.list_setups(day.id), key=lambda value: value.order):
                setup_label = setup.name + (f" / {setup.setup_code}" if setup.setup_code else "")
                setup_item = f"setup:{setup.id}"
                self.tree.insert(day_item, "end", iid=setup_item, text=setup_label, open=True)
                laps = self.services.list_laps(setup.id)
                best_by_event = self.services.best_laps_by_event(setup.id)
                pending_order = None
                if (
                    self._current_kind == "lap"
                    and self._pending_lap_order is not None
                    and isinstance(self._record, Lap)
                    and self._record.setup_id == setup.id
                ):
                    pending_order = self._pending_lap_order
                if pending_order is None:
                    ordered_laps = sorted(laps, key=lambda value: value.sequence)
                else:
                    positions = {lap_id: index for index, lap_id in enumerate(pending_order)}
                    ordered_laps = sorted(
                        laps,
                        key=lambda value: (positions.get(value.id, len(positions) + value.sequence), value.sequence),
                    )
                for lap in ordered_laps:
                    event = self.services.get_event_layout(lap.event_layout_id)
                    event_label = (
                        f"{event.event_name + ' · ' if event.event_name else ''}{event.track_name} / {event.layout_name}"
                        if event is not None else "Unknown event/layout"
                    )
                    label = f"{format_lap_time(lap.time_ms)} · {event_label}"
                    if lap.status == "invalid":
                        label += " · INVALID"
                    best = best_by_event.get(lap.event_layout_id)
                    if lap.status == "valid" and best is not None and best.id == lap.id:
                        label += " · BEST"
                    self.tree.insert(setup_item, "end", iid=f"lap:{lap.id}", text=label)
        if desired and self.tree.exists(desired):
            self.tree.selection_set(desired)
            self.tree.see(desired)
            self._open_tree_ancestors(desired)
        self._suppress_tree_selection = False
        self._refresh_event_choices()
        self._update_action_states()

    def _open_tree_ancestors(self, item: str) -> None:
        ancestor = self.tree.parent(item)
        while ancestor:
            self.tree.item(ancestor, open=True)
            ancestor = self.tree.parent(ancestor)

    def _on_tree_select(self, _event=None) -> None:
        if self._suppress_tree_selection or self.tree is None:
            return
        selected = self.tree.selection()
        if not selected:
            return
        item = selected[0]
        target_kind, separator, identifier = item.partition(":")
        if not separator or target_kind not in ("day", "setup", "lap"):
            return
        if (
            not self._is_new and self._current_kind == target_kind
            and self._record is not None and self._record.id == identifier
        ):
            return
        if not self._resolve_unsaved():
            self._restore_tree_selection()
            return
        getter = {"day": self.services.get_day, "setup": self.services.get_setup, "lap": self.services.get_lap}[target_kind]
        record = getter(identifier)
        if record is None:
            self._set_status("This record no longer exists.")
            self.refresh_tree()
            return
        self._show_editor(target_kind, record, is_new=False)
        self._restore_tree_selection()

    def _restore_tree_selection(self) -> None:
        if self.tree is None:
            return
        previous = self._tree_item_for_current()
        self._suppress_tree_selection = True
        if previous and self.tree.exists(previous):
            self.tree.selection_set(previous)
            self.tree.see(previous)
            self._open_tree_ancestors(previous)
        else:
            self.tree.selection_remove(self.tree.selection())
        self._suppress_tree_selection = False

    def _has_pending_lap_order(self) -> bool:
        return (
            self._current_kind == "lap"
            and self._pending_lap_order is not None
            and self._pending_lap_order != self._lap_order_original
        )

    def _is_dirty(self) -> bool:
        if self._staged:
            return True
        if self._has_pending_lap_order():
            return True
        values = self._capture_values()
        if self._is_new:
            return any(value.strip() for value in values.values())
        return values != self._loaded_values

    def _update_dirty_indicator(self) -> None:
        if self._has_pending_lap_order():
            dirty_text = "Unsaved changes · lap order pending"
        else:
            dirty_text = "Unsaved changes" if self._is_dirty() else ""
        self._set_dirty_text(dirty_text)
        self._update_action_states()

    def _discard_staged(self) -> None:
        if not self._staged:
            return
        failures = self.services.discard_staged(tuple(self._staged))
        self._staged.clear()
        self._refresh_attachment_list()
        if failures:
            message = "Could not remove staged file copies:\n" + "\n".join(map(str, failures))
            self._set_status(message)
            messagebox.showwarning("Staged files need cleanup", message, parent=self.root)

    def _discard_pending_lap_order(self) -> None:
        if not self._has_pending_lap_order():
            return
        self._pending_lap_order = list(self._lap_order_original)
        self.refresh_tree()
        self._update_dirty_indicator()

    def _resolve_unsaved(self) -> bool:
        if not self._is_dirty():
            return True
        answer = messagebox.askyesnocancel(
            "Unsaved changes",
            "Save your changes before continuing?\nChoose No to discard them.",
            parent=self.root,
        )
        if answer is None:
            return False
        if answer:
            return self.save_current()
        self._discard_staged()
        self._discard_pending_lap_order()
        return True

    def _selection_item(self) -> str | None:
        if self.tree is None:
            return None
        selection = self.tree.selection()
        return selection[0] if selection else None


    def _selected_setup_id(self) -> str | None:
        if self._current_kind == "setup":
            return self._record.id if isinstance(self._record, Setup) and not self._is_new else None
        if self._current_kind == "lap":
            return self._record.setup_id if isinstance(self._record, Lap) else self._context_id
        item = self._selection_item()
        if item:
            kind, _, identifier = item.partition(":")
            if kind == "setup":
                return identifier
            if kind == "lap":
                lap = self.services.get_lap(identifier)
                return lap.setup_id if lap is not None else None
        return None

    def _selected_day_id(self) -> str | None:
        if self._current_kind == "day":
            return self._record.id if isinstance(self._record, TestDay) and not self._is_new else None
        if self._current_kind == "setup":
            return self._record.test_day_id if isinstance(self._record, Setup) else self._context_id
        setup_id = self._selected_setup_id()
        if setup_id:
            setup = self.services.get_setup(setup_id)
            return setup.test_day_id if setup is not None else None
        item = self._selection_item()
        if item:
            kind, _, identifier = item.partition(":")
            if kind == "day":
                return identifier
            if kind == "setup":
                setup = self.services.get_setup(identifier)
                return setup.test_day_id if setup is not None else None
            if kind == "lap":
                lap = self.services.get_lap(identifier)
                setup = self.services.get_setup(lap.setup_id) if lap is not None else None
                return setup.test_day_id if setup is not None else None
        return None

    def _update_action_states(self) -> None:
        if not self.actions:
            return
        day_id = self._selected_day_id()
        setup_id = self._selected_setup_id()
        record_saved = self._record is not None and not self._is_new
        has_lap = self._current_kind == "lap" and record_saved
        self.actions["add_setup"].state(["!disabled"] if day_id else ["disabled"])
        self.actions["add_lap"].state(["!disabled"] if setup_id else ["disabled"])
        duplicate_ok = bool(setup_id and self.services.get_setup(setup_id))
        self.actions["duplicate_setup"].state(["!disabled"] if duplicate_ok else ["disabled"])
        self.actions["move_lap_up"].state(["!disabled"] if has_lap else ["disabled"])
        self.actions["move_lap_down"].state(["!disabled"] if has_lap else ["disabled"])
        self.actions["save"].state(["!disabled"] if self._current_kind else ["disabled"])
        self.actions["delete"].state(["!disabled"] if record_saved else ["disabled"])
        self.actions["attach_file"].state(["!disabled"] if self._current_kind else ["disabled"])
        self.actions["open_attachment"].state(["!disabled"] if self._attachment_refs else ["disabled"])
        self.actions["export_csv"].state(["!disabled"] if setup_id else ["disabled"])

    def add_day(self) -> None:
        if self._resolve_unsaved():
            self._show_editor("day", None, is_new=True)

    def add_setup(self) -> None:
        if not self._resolve_unsaved():
            return
        day_id = self._selected_day_id()
        if day_id is None:
            self._set_status("Select a saved test day before adding a setup.")
            return
        self._show_editor("setup", None, is_new=True, context_id=day_id)

    def add_lap(self) -> None:
        if not self._resolve_unsaved():
            return
        setup_id = self._selected_setup_id()
        if setup_id is None:
            self._set_status("Select a saved setup before adding a lap.")
            return
        self._show_editor("lap", None, is_new=True, context_id=setup_id)

    def duplicate_setup(self) -> None:
        if not self._resolve_unsaved():
            return
        setup_id = self._selected_setup_id()
        if setup_id is None:
            self._set_status("Select a saved setup to duplicate.")
            return
        try:
            duplicate = self.services.duplicate_setup(setup_id)
        except Exception as error:
            self._set_status(f"Could not duplicate setup: {error}")
            messagebox.showerror("Duplicate setup", str(error), parent=self.root)
            return
        self.refresh_tree(select_item=f"setup:{duplicate.id}")
        self._show_editor("setup", duplicate, is_new=False)
        self._set_status("Setup duplicated. Review its details and save any changes.")

    def _new_record_values(self, kind: str):
        identifier = new_id()
        timestamp = utc_now_iso()
        values = self._capture_values()
        if kind == "day":
            existing = self._record if isinstance(self._record, TestDay) else None
            return TestDay(
                id=existing.id if existing else identifier,
                date=values["date"],
                location=values["location"],
                weather=self._optional(values["weather"]),
                notes=self._optional(values["notes"]),
                created_at=existing.created_at if existing else timestamp,
            )
        if kind == "setup":
            existing_record = self._record if isinstance(self._record, Setup) else None
            day_id = existing_record.test_day_id if existing_record else self._context_id
            if day_id is None:
                raise ValidationError("test_day_id", "Select a test day for this setup.")
            existing = self.services.list_setups(day_id)
            return Setup(
                id=existing_record.id if existing_record else identifier,
                test_day_id=day_id,
                name=values["name"],
                setup_code=self._optional(values["setup_code"]),
                settings_text=values["settings_text"],
                notes=self._optional(values["notes"]),
                order=existing_record.order if existing_record else len(existing) + 1,
                created_at=existing_record.created_at if existing_record else timestamp,
            )
        if kind == "lap":
            existing_record = self._record if isinstance(self._record, Lap) else None
            setup_id = existing_record.setup_id if existing_record else self._context_id
            if setup_id is None:
                raise ValidationError("setup_id", "Select a setup for this lap.")
            event_id = self._event_id_by_label.get(values["event_layout_id"])
            if not event_id:
                raise ValidationError("event_layout_id", "Choose an active event/layout.")
            existing = self.services.list_laps(setup_id)
            return Lap(
                id=existing_record.id if existing_record else identifier,
                setup_id=setup_id,
                event_layout_id=event_id,
                sequence=existing_record.sequence if existing_record else len(existing) + 1,
                time_ms=parse_lap_time(values["lap_time"]),
                status=values["status"],
                driver=self._optional(values["driver"]),
                time_of_day=self._optional(values["time_of_day"]),
                notes=self._optional(values["notes"]),
                created_at=existing_record.created_at if existing_record else timestamp,
            )
        raise ValueError(f"Unknown record kind: {kind}")

    def save_current(self) -> bool:
        if self._current_kind is None:
            return False
        self._clear_field_errors()
        try:
            record = self._new_record_values(self._current_kind)
        except ValidationError as error:
            self._show_validation_error(error)
            return False
        except Exception as error:
            self._set_status(f"Could not read form values: {error}")
            return False
        try:
            if self._current_kind == "lap":
                lap_order = self._pending_lap_order
                if lap_order == self._lap_order_original:
                    lap_order = None
                saved = self.services.save_lap(
                    record,
                    tuple(self._staged),
                    lap_order=tuple(lap_order) if lap_order is not None else None,
                )
            else:
                saved = getattr(self.services, f"save_{self._current_kind}")(record, tuple(self._staged))
        except ValidationError as error:
            had_staged = bool(self._staged)
            self._staged.clear()
            self._refresh_attachment_list()
            self._show_validation_error(error)
            if had_staged:
                self._set_status(f"{error.message} Staged file copies were removed; reattach them before retrying.")
            return False
        except Exception as error:
            had_staged = bool(self._staged)
            self._staged.clear()
            self._refresh_attachment_list()
            message = f"Save failed: {error}"
            if had_staged:
                message += " Staged file copies were removed; reattach them before retrying."
            self._set_status(message)
            messagebox.showerror(f"Save {self._current_kind}", message, parent=self.root)
            return False

        kind = self._current_kind
        self._staged.clear()
        self._is_new = False
        self._record = saved
        self._context_id = None
        self.refresh_tree(select_item=f"{kind}:{saved.id}")
        self._show_editor(kind, saved, is_new=False)
        self._set_status(f"{kind.title()} saved.")
        return True


    def _refresh_attachment_list(self) -> None:
        if self.attachment_list is None:
            return
        self.attachment_list.delete(0, "end")
        self._attachment_refs = []
        if self._record is not None and not self._is_new and self._current_kind is not None:
            for attachment in self.services.list_attachments(self._current_kind, self._record.id):
                self._attachment_refs.append(("saved", attachment))
                self.attachment_list.insert("end", attachment.original_name)
        for staged in self._staged:
            self._attachment_refs.append(("staged", staged))
            self.attachment_list.insert("end", f"{staged.original_name} (staged)")
        self._update_action_states()

    def attach_file(self) -> None:
        if self._current_kind is None:
            self._set_status("Select or add a test day, setup, or lap before attaching a file.")
            return
        source = filedialog.askopenfilename(parent=self.root, title="Choose an attachment")
        if not source:
            return
        try:
            staged = self.services.stage_attachment(Path(source), "file")
        except Exception as error:
            self._set_status(f"Could not stage attachment: {error}")
            messagebox.showerror("Attach file", str(error), parent=self.root)
            return
        self._staged.append(staged)
        self._refresh_attachment_list()
        self._update_dirty_indicator()
        self._set_status("File copied into the data folder. Save this record to keep it.")

    def open_data_folder(self) -> None:
        try:
            opener = getattr(os, "startfile", None)
            if opener is None:
                raise OSError("Opening folders is available on Windows.")
            opener(str(self.services.paths.root))
        except Exception as error:
            self._set_status(f"Could not open data folder: {error}")
            messagebox.showerror("Open data folder", str(error), parent=self.root)

    def open_attachment(self) -> None:
        selection = self.attachment_list.curselection() if self.attachment_list is not None else ()
        if not selection:
            self._set_status("Select a saved attachment to open.")
            return
        ref_kind, attachment = self._attachment_refs[selection[0]]
        if ref_kind == "staged":
            self._set_status("Save the record before opening a newly attached file.")
            return
        try:
            path = self.services.resolve_attachment(attachment)
            opener = getattr(os, "startfile", None)
            if opener is None:
                raise OSError("Opening files is available on Windows.")
            opener(str(path))
        except Exception as error:
            self._set_status(f"Could not open attachment: {error}")
            messagebox.showerror("Open attachment", str(error), parent=self.root)

    def manage_events(self) -> None:
        if not self._resolve_unsaved():
            return
        if self.event_manager is not None:
            try:
                if self.event_manager.winfo_exists():
                    self.event_manager.lift()
                    return
            except tk.TclError:
                pass
        self.event_manager = EventLayoutManager(self.root, self.services, on_change=self.refresh_tree)

    def _current_saved_record(self):
        if self._record is not None and not self._is_new:
            return self._current_kind, self._record
        return None, None

    def _item_from_selection(self):
        item = self._selection_item()
        if not item:
            return None, None
        kind, separator, identifier = item.partition(":")
        if not separator:
            return None, None
        getter = {
            "day": self.services.get_day,
            "setup": self.services.get_setup,
            "lap": self.services.get_lap,
        }.get(kind)
        return (kind, getter(identifier)) if getter else (None, None)

    def delete_current(self) -> None:
        if not self._resolve_unsaved():
            return
        kind, record = self._current_saved_record()
        if record is None:
            kind, record = self._item_from_selection()
        if record is None or kind not in ("day", "setup", "lap"):
            self._set_status("Select a saved test day, setup, or lap to delete.")
            return
        if kind == "day":
            setups = self.services.list_setups(record.id)
            lap_count = sum(len(self.services.list_laps(setup.id)) for setup in setups)
            prompt = (
                f"Delete test day {record.date} at {record.location}? "
                f"This also removes {len(setups)} setup(s), {lap_count} lap(s), and their attachment records."
            )
            title = "Delete test day"
        elif kind == "setup":
            laps = self.services.list_laps(record.id)
            prompt = (
                f"Delete setup {record.name}? This also removes {len(laps)} lap(s) "
                "and the setup's and laps' attachment records."
            )
            title = "Delete setup"
        else:
            prompt = f"Delete lap {format_lap_time(record.time_ms)} and its attachment records?"
            title = "Delete lap"
        if not messagebox.askyesno(title, prompt, parent=self.root):
            return
        try:
            failures = getattr(self.services, f"delete_{kind}")(record.id)
        except Exception as error:
            self._set_status(f"Could not delete {kind}: {error}")
            messagebox.showerror(title, str(error), parent=self.root)
            return
        self._discard_staged()
        self._show_editor(None, None, is_new=True)
        self.refresh_tree()
        if failures:
            messagebox.showwarning(
                "Copied files need cleanup",
                "The records were deleted, but these copied files could not be removed:\n"
                + "\n".join(map(str, failures)),
                parent=self.root,
            )
        else:
            self._set_status(f"{kind.title()} deleted.")

    def move_lap(self, direction: int) -> None:
        kind, record = self._current_saved_record()
        if kind != "lap" or not isinstance(record, Lap):
            kind, record = self._item_from_selection()
        if kind != "lap" or not isinstance(record, Lap):
            self._set_status("Select a saved lap to change its order.")
            return
        if (
            self._current_kind == "lap"
            and self._pending_lap_order is not None
            and isinstance(self._record, Lap)
            and self._record.setup_id == record.setup_id
        ):
            lap_ids = list(self._pending_lap_order)
        else:
            lap_ids = [lap.id for lap in self.services.list_laps(record.setup_id)]
        try:
            current_index = lap_ids.index(record.id)
        except ValueError:
            return
        destination = current_index + direction
        if destination < 0 or destination >= len(lap_ids):
            return
        lap_ids[current_index], lap_ids[destination] = lap_ids[destination], lap_ids[current_index]
        if self._current_kind == "lap" and isinstance(self._record, Lap):
            self._pending_lap_order = lap_ids
        else:
            self._lap_order_original = [lap.id for lap in self.services.list_laps(record.setup_id)]
            self._pending_lap_order = lap_ids
        self.refresh_tree(select_item=f"lap:{record.id}")
        self._update_dirty_indicator()
        if self._pending_lap_order != self._lap_order_original:
            self._set_status("Lap order pending. Save to apply or choose Discard when leaving this lap.")
        else:
            self._set_status("Lap order restored to its saved order.")

    def create_backup(self) -> None:
        if not self._resolve_unsaved():
            return
        try:
            destination = self.services.create_backup()
        except Exception as error:
            self._set_status(f"Backup failed: {error}")
            messagebox.showerror("Create backup", str(error), parent=self.root)
            return
        self._set_status(f"Backup created: {destination}")

    def export_csv(self) -> None:
        if not self._resolve_unsaved():
            return
        setup_id = self._selected_setup_id()
        if setup_id is None:
            self._set_status("Select a setup or one of its laps before exporting CSV.")
            return
        destination = filedialog.asksaveasfilename(
            parent=self.root,
            title="Export setup laps",
            defaultextension=".csv",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            initialfile="test-laps.csv",
        )
        if not destination:
            return
        try:
            written = self.services.export_setup(setup_id, Path(destination))
        except Exception as error:
            self._set_status(f"CSV export failed: {error}")
            messagebox.showerror("Export CSV", str(error), parent=self.root)
            return
        self._set_status(f"CSV exported: {written}")

    def close_request(self) -> bool:
        if not self._resolve_unsaved():
            return False
        if self.event_manager is not None:
            try:
                if self.event_manager.winfo_exists() and not self.event_manager.close_request():
                    return False
            except tk.TclError:
                pass
        if not self._closed:
            self._closed = True
            self.services.close()
        self.root.destroy()
        return True
