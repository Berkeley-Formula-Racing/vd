"""Event and track-layout library window for the offline test log."""

from __future__ import annotations

import math
import os
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from .models import EventLayout, StagedAttachment, new_id, utc_now_iso
from .validation import ValidationError


class EventLayoutManager(tk.Toplevel):
    """Create, duplicate, archive, and edit reusable event/layout records."""

    def __init__(self, parent, services, on_change=None):
        super().__init__(parent)
        self.parent = parent
        self.services = services
        self.on_change = on_change
        self.title("Manage Events and Layouts")
        self.geometry("980x560")
        self.minsize(800, 450)

        self.tree = None
        self.fields: dict[str, tk.Widget] = {}
        self.actions: dict[str, ttk.Button] = {}
        self.attachment_list: tk.Listbox | None = None
        self.status_label: ttk.Label | None = None
        self._current_item: str | None = None
        self._record: EventLayout | None = None
        self._is_new = False
        self._loaded_values: dict[str, str] = {}
        self._staged: list[StagedAttachment] = []
        self._attachments = []
        self._suppress_selection = False
        self._errors: dict[str, ttk.Label] = {}
        self.protocol("WM_DELETE_WINDOW", self.close_request)
        self._build_widgets()
        self.refresh_tree()
        self._show_blank()

    def _build_widgets(self) -> None:
        toolbar = ttk.Frame(self, padding=(8, 8, 8, 4))
        toolbar.pack(fill="x")
        button_specs = (
            ("new_event", "New Event/Layout", self.new_event),
            ("duplicate", "Duplicate", self.duplicate_event),
            ("save", "Save", self.save_event),
            ("archive", "Archive / Restore", self.toggle_archive),
            ("delete", "Delete", self.delete_event),
            ("attach_map", "Attach Map", self.attach_map),
            ("open_map", "Open Map", self.open_map),
        )
        for column, (name, label, callback) in enumerate(button_specs):
            button = ttk.Button(toolbar, text=label, command=callback)
            button.grid(row=0, column=column, padx=(0, 5), sticky="ew")
            toolbar.columnconfigure(column, weight=1)
            self.actions[name] = button

        body = ttk.Panedwindow(self, orient="horizontal")
        body.pack(fill="both", expand=True, padx=8, pady=(4, 8))
        tree_frame = ttk.Frame(body, padding=(0, 0, 8, 0))
        editor_frame = ttk.Frame(body)
        body.add(tree_frame, weight=1)
        body.add(editor_frame, weight=2)

        self.tree = ttk.Treeview(tree_frame, show="tree", selectmode="browse")
        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.tree.bind("<<TreeviewSelect>>", self._on_tree_select)

        editor_frame.columnconfigure(1, weight=1)
        self._editor_frame = editor_frame
        form = ttk.Frame(editor_frame, padding=(8, 4))
        form.grid(row=0, column=0, sticky="nsew")
        form.columnconfigure(1, weight=1)
        form.columnconfigure(2, weight=0)
        self._form_frame = form
        self._add_entry(form, "track_name", "Track / venue", 0)
        self._add_entry(form, "layout_name", "Layout", 1)
        self._add_entry(form, "event_name", "Event name", 2)
        self._add_entry(form, "event_type", "Event type", 3)
        self._add_entry(form, "length_m", "Length (m)", 4)
        ttk.Label(form, text="Notes").grid(row=5, column=0, sticky="nw", padx=(0, 8), pady=4)
        notes = tk.Text(form, height=4, width=38, wrap="word")
        notes.grid(row=5, column=1, sticky="ew", pady=4)
        self.fields["notes"] = notes
        self._errors["notes"] = ttk.Label(form, text="", foreground="#a02020")
        self._errors["notes"].grid(row=6, column=1, sticky="w")

        attachments_frame = ttk.LabelFrame(editor_frame, text="Map files", padding=6)
        attachments_frame.grid(row=1, column=0, sticky="nsew", padx=8, pady=(2, 6))
        self.attachment_list = tk.Listbox(attachments_frame, height=5, exportselection=False)
        self.attachment_list.pack(fill="both", expand=True)
        self.status_label = ttk.Label(editor_frame, text="", foreground="#7a2525", wraplength=520)
        self.status_label.grid(row=2, column=0, sticky="ew", padx=8, pady=(0, 6))
        editor_frame.rowconfigure(1, weight=1)

    def _add_entry(self, parent, key: str, label: str, row: int) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=4)
        widget = ttk.Entry(parent)
        widget.grid(row=row, column=1, sticky="ew", pady=4)
        self.fields[key] = widget
        self._errors[key] = ttk.Label(parent, text="", foreground="#a02020")
        self._errors[key].grid(row=row, column=2, sticky="w")

    def _capture_values(self) -> dict[str, str]:
        values = {}
        for key, widget in self.fields.items():
            if isinstance(widget, tk.Text):
                values[key] = widget.get("1.0", "end-1c")
            else:
                values[key] = widget.get()
        return values

    def _set_values(self, values: dict[str, str]) -> None:
        for key, widget in self.fields.items():
            value = values.get(key, "")
            if isinstance(widget, tk.Text):
                widget.configure(state="normal")
                widget.delete("1.0", "end")
                widget.insert("1.0", value)
            else:
                widget.configure(state="normal")
                widget.delete(0, "end")
                widget.insert(0, value)
        self._loaded_values = self._capture_values()
        self._clear_errors()

    def _display_values(self, record: EventLayout) -> dict[str, str]:
        return {
            "track_name": record.track_name,
            "layout_name": record.layout_name,
            "event_name": record.event_name or "",
            "event_type": record.event_type or "",
            "length_m": "" if record.length_m is None else str(record.length_m),
            "notes": record.notes or "",
        }

    def _clear_errors(self) -> None:
        for label in self._errors.values():
            label.configure(text="")

    def _show_error(self, field: str, message: str) -> None:
        if field in self._errors:
            self._errors[field].configure(text=message)
        self._set_status(message)

    def _set_status(self, text: str) -> None:
        if self.status_label is not None:
            self.status_label.configure(text=text)

    def _set_editable(self, editable: bool) -> None:
        state = "normal" if editable else "disabled"
        for widget in self.fields.values():
            widget.configure(state=state)
        self.actions["save"].state(["!disabled"] if editable else ["disabled"])
        self.actions["attach_map"].state(["!disabled"] if editable else ["disabled"])
        self.actions["archive"].state(["!disabled"] if self._record is not None else ["disabled"])
        self.actions["delete"].state(
            ["!disabled"] if self._record is not None and not self.services.event_is_referenced(self._record.id)
            else ["disabled"]
        )
        self.actions["duplicate"].state(["!disabled"] if self._record is not None else ["disabled"])
        self.actions["open_map"].state(["!disabled"] if self._attachments else ["disabled"])

    def refresh_tree(self, select_id: str | None = None) -> None:
        if self.tree is None:
            return
        previous_id = select_id or (self._record.id if self._record is not None else None)
        self._suppress_selection = True
        for item in self.tree.get_children(""):
            self.tree.delete(item)
        records = self.services.list_event_layouts(include_archived=True)
        records.sort(key=lambda item: (item.track_name.casefold(), item.layout_name.casefold(), item.event_name or ""))
        for record in records:
            label = f"{record.event_name + ' · ' if record.event_name else ''}{record.track_name} / {record.layout_name}"
            if record.archived:
                label += " [Archived]"
            self.tree.insert("", "end", iid=f"event:{record.id}", text=label)
        if previous_id and self.tree.exists(f"event:{previous_id}"):
            self.tree.selection_set(f"event:{previous_id}")
            self.tree.see(f"event:{previous_id}")
        self._suppress_selection = False

    def _show_blank(self) -> None:
        self._current_item = None
        self._record = None
        self._is_new = True
        self._attachments = []
        self._staged.clear()
        self._set_values({key: "" for key in self.fields})
        self._refresh_attachment_list()
        self._set_status("")
        self._set_editable(True)
        self.actions["delete"].state(["disabled"])
        self.actions["archive"].state(["disabled"])
        self.actions["duplicate"].state(["disabled"])
        self.actions["open_map"].state(["disabled"])

    def _load_event(self, identifier: str) -> None:
        record = self.services.get_event_layout(identifier)
        if record is None:
            self._show_blank()
            self._set_status("This event/layout no longer exists.")
            return
        self._current_item = f"event:{record.id}"
        self._record = record
        self._is_new = False
        self._staged.clear()
        self._attachments = self.services.list_attachments("event_layout", record.id)
        self._set_values(self._display_values(record))
        self._refresh_attachment_list()
        self._set_status("Archived records remain on historical laps." if record.archived else "")
        editable = not self.services.event_is_referenced(record.id)
        self._set_editable(editable)

    def _on_tree_select(self, _event=None) -> None:
        if self._suppress_selection:
            return
        selection = self.tree.selection()
        if not selection:
            return
        item = selection[0]
        if item == self._current_item:
            return
        if not self._resolve_unsaved():
            self._restore_selection()
            return
        self._load_event(item.partition(":")[2])

    def _restore_selection(self) -> None:
        self._suppress_selection = True
        if self._current_item and self.tree.exists(self._current_item):
            self.tree.selection_set(self._current_item)
        else:
            self.tree.selection_remove(self.tree.selection())
        self._suppress_selection = False

    def _is_dirty(self) -> bool:
        if self._staged:
            return True
        values = self._capture_values()
        if self._is_new:
            return any(value.strip() for value in values.values())
        return values != self._loaded_values

    def _resolve_unsaved(self) -> bool:
        if not self._is_dirty():
            return True
        answer = messagebox.askyesnocancel(
            "Unsaved event/layout changes",
            "Save your event/layout changes before continuing?\nChoose No to discard them.",
            parent=self,
        )
        if answer is None:
            return False
        if answer:
            return self.save_event()
        self._discard_staged()
        return True

    def _discard_staged(self) -> None:
        if self._staged:
            failures = self.services.discard_staged(tuple(self._staged))
            self._staged.clear()
            self._refresh_attachment_list()
            if failures:
                self._set_status("Could not remove staged map copies:\n" + "\n".join(map(str, failures)))

    def _clear_tree_selection(self) -> None:
        self._suppress_selection = True
        self.tree.selection_remove(self.tree.selection())
        self._suppress_selection = False

    def new_event(self) -> None:
        if not self._resolve_unsaved():
            self._restore_selection()
            return
        self._clear_tree_selection()
        self._show_blank()

    def duplicate_event(self) -> None:
        if not self._resolve_unsaved() or self._record is None:
            return
        try:
            duplicate = self.services.duplicate_event_layout(self._record.id)
        except Exception as error:
            self._set_status(f"Could not duplicate event/layout: {error}")
            messagebox.showerror("Duplicate event/layout", str(error), parent=self)
            return
        self.refresh_tree(select_id=duplicate.id)
        self._load_event(duplicate.id)
        self._set_status("Event/layout duplicated. Update its details and save.")

    def _optional(self, value: str) -> str | None:
        value = value.strip()
        return value or None

    def _make_record(self) -> EventLayout:
        values = self._capture_values()
        length_text = values["length_m"].strip()
        if length_text:
            try:
                length = float(length_text)
            except ValueError as error:
                raise ValidationError("length_m", "Enter a positive length in metres.") from error
            if not math.isfinite(length) or length <= 0:
                raise ValidationError("length_m", "Length must be a finite positive number of metres.")
        else:
            length = None
        return EventLayout(
            id=self._record.id if self._record is not None else new_id(),
            track_name=values["track_name"],
            layout_name=values["layout_name"],
            event_name=self._optional(values["event_name"]),
            event_type=self._optional(values["event_type"]),
            length_m=length,
            notes=self._optional(values["notes"]),
            archived=self._record.archived if self._record is not None else False,
            created_at=self._record.created_at if self._record is not None else utc_now_iso(),
        )

    def save_event(self) -> bool:
        self._clear_errors()
        if self._record is not None and self.services.event_is_referenced(self._record.id):
            self._set_status("Referenced event/layout details are frozen. Duplicate it to make changes.")
            return False
        try:
            record = self._make_record()
        except ValidationError as error:
            self._show_error(error.field, error.message)
            return False
        except Exception as error:
            self._set_status(f"Could not read event/layout fields: {error}")
            return False
        try:
            saved = self.services.save_event_layout(record, tuple(self._staged))
        except Exception as error:
            had_staged = bool(self._staged)
            self._staged.clear()
            self._refresh_attachment_list()
            message = f"Save failed: {error}"
            if had_staged:
                message += " Staged map copies were removed; reattach them before retrying."
            self._set_status(message)
            messagebox.showerror("Save event/layout", message, parent=self)
            return False
        self._staged.clear()
        self._record = saved
        self._is_new = False
        self._current_item = f"event:{saved.id}"
        self._attachments = self.services.list_attachments("event_layout", saved.id)
        self._set_values(self._display_values(saved))
        self._refresh_attachment_list()
        self.refresh_tree(select_id=saved.id)
        self._set_editable(not self.services.event_is_referenced(saved.id))
        self._set_status("Event/layout saved.")
        if self.on_change is not None:
            self.on_change()
        return True

    def toggle_archive(self) -> None:
        if not self._resolve_unsaved() or self._record is None:
            return
        try:
            saved = self.services.archive_event_layout(self._record.id, not self._record.archived)
        except Exception as error:
            self._set_status(f"Could not archive event/layout: {error}")
            messagebox.showerror("Archive event/layout", str(error), parent=self)
            return
        self._record = saved
        self._is_new = False
        self._current_item = f"event:{saved.id}"
        self._set_values(self._display_values(saved))
        self.refresh_tree(select_id=saved.id)
        self._set_editable(not self.services.event_is_referenced(saved.id))
        self._set_status("Event/layout archived." if saved.archived else "Event/layout restored.")
        if self.on_change is not None:
            self.on_change()

    def delete_event(self) -> None:
        if self._record is None:
            return
        if self.services.event_is_referenced(self._record.id):
            self._set_status("This event/layout is used by a lap and cannot be deleted.")
            return
        label = f"{self._record.event_name or self._record.track_name} / {self._record.layout_name}"
        if not messagebox.askyesno(
            "Delete event/layout",
            f"Delete event/layout “{label}” and its map attachments?",
            parent=self,
        ):
            return
        try:
            failures = self.services.delete_event_layout(self._record.id)
        except Exception as error:
            self._set_status(f"Could not delete event/layout: {error}")
            messagebox.showerror("Delete event/layout", str(error), parent=self)
            return
        self._show_blank()
        self.refresh_tree()
        if failures:
            messagebox.showwarning(
                "Map files need cleanup",
                "The event/layout was deleted, but these copied files could not be removed:\n"
                + "\n".join(map(str, failures)),
                parent=self,
            )
        else:
            self._set_status("Event/layout deleted.")
        if self.on_change is not None:
            self.on_change()

    def attach_map(self) -> None:
        if self._record is not None and self.services.event_is_referenced(self._record.id):
            self._set_status("Maps on referenced event/layouts are frozen. Duplicate the record to add a map.")
            return
        source = filedialog.askopenfilename(parent=self, title="Choose a map file")
        if not source:
            return
        try:
            staged = self.services.stage_attachment(Path(source), "map")
        except Exception as error:
            self._set_status(f"Could not stage map file: {error}")
            messagebox.showerror("Attach map", str(error), parent=self)
            return
        self._staged.append(staged)
        self._refresh_attachment_list()
        self._set_status("Map copied into the data folder. Save the event/layout to keep it.")

    def _refresh_attachment_list(self) -> None:
        if self.attachment_list is None:
            return
        self.attachment_list.delete(0, "end")
        for attachment in self._attachments:
            self.attachment_list.insert("end", attachment.original_name)
        for attachment in self._staged:
            self.attachment_list.insert("end", f"{attachment.original_name} (staged)")

    def open_map(self) -> None:
        selection = self.attachment_list.curselection() if self.attachment_list is not None else ()
        if not selection:
            self._set_status("Select a saved map to open.")
            return
        index = selection[0]
        if index >= len(self._attachments):
            self._set_status("Save the event/layout before opening a newly attached map.")
            return
        attachment = self._attachments[index]
        try:
            path = self.services.resolve_attachment(attachment)
            opener = getattr(os, "startfile", None)
            if opener is None:
                raise OSError("Opening files is available on Windows.")
            opener(str(path))
        except Exception as error:
            self._set_status(f"Could not open map: {error}")
            messagebox.showerror("Open map", str(error), parent=self)

    def _refresh_buttons(self) -> None:
        has_record = self._record is not None
        if not has_record:
            self.actions["duplicate"].state(["disabled"])
            self.actions["archive"].state(["disabled"])
            self.actions["delete"].state(["disabled"])
            self.actions["open_map"].state(["disabled"])
            return
        referenced = self.services.event_is_referenced(self._record.id)
        self.actions["duplicate"].state(["!disabled"])
        self.actions["archive"].state(["!disabled"])
        self.actions["delete"].state(["disabled"] if referenced else ["!disabled"])
        self.actions["open_map"].state(["!disabled"] if self._attachments else ["disabled"])

    def close_request(self) -> bool:
        if not self._resolve_unsaved():
            return False
        self.destroy()
        return True
