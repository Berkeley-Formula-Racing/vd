"""Callback tests for the Tkinter test log interface."""

import importlib
import tempfile
import tkinter as tk
from tkinter import ttk
import unittest
from pathlib import Path
from unittest.mock import patch

from vd_test_log.models import DataPaths, EventLayout, Lap, Setup, TestDay, new_id, utc_now_iso
from vd_test_log.services import TestLogServices


class TestLogUiTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root_path = Path(self.temporary_directory.name)
        data_root = self.root_path / "data"
        self.paths = DataPaths(
            root=data_root,
            database=data_root / "test_log.sqlite3",
            attachments=data_root / "attachments",
            lock_file=data_root / ".test_log.lock",
            backup_root=data_root / "backups",
        )
        self.services = TestLogServices.open(self.paths)
        self.root = tk.Tk()
        self.root.withdraw()
        self.ui_module = importlib.import_module("vd_test_log.ui")
        self.window = None
        self.managers = []

    def tearDown(self):
        for manager in self.managers:
            try:
                if manager.winfo_exists():
                    manager.destroy()
            except tk.TclError:
                pass
        try:
            self.root.destroy()
        except tk.TclError:
            pass
        self.services.close()
        self.temporary_directory.cleanup()

    def require_window(self):
        window_type = getattr(self.ui_module, "TestLogWindow", None)
        self.assertIsNotNone(window_type, "vd_test_log.ui must expose TestLogWindow")
        self.window = window_type(self.root, self.services)
        self.window.build()
        self.root.update()
        return self.window

    def make_day(self, location="Test Site"):
        return TestDay(
            id=new_id(),
            date="2026-10-05",
            location=location,
            weather=None,
            notes=None,
            created_at=utc_now_iso(),
        )

    def make_setup(self, day_id, name="Baseline", order=1):
        return Setup(
            id=new_id(),
            test_day_id=day_id,
            name=name,
            setup_code="A01",
            settings_text="Cold pressures: 20 psi",
            notes="Keep the front settings",
            order=order,
            created_at=utc_now_iso(),
        )

    def make_event(self, track="Test Site", layout="Short Loop", event_name="Handling Day"):
        return EventLayout(
            id=new_id(),
            track_name=track,
            layout_name=layout,
            event_name=event_name,
            event_type="test",
            length_m=1234.5,
            notes="Dry surface",
            archived=False,
            created_at=utc_now_iso(),
        )

    def make_lap(self, setup_id, event_id, sequence=1, time_ms=42318):
        return Lap(
            id=new_id(),
            setup_id=setup_id,
            event_layout_id=event_id,
            sequence=sequence,
            time_ms=time_ms,
            status="valid",
            driver="Driver",
            time_of_day=None,
            notes=None,
            created_at=utc_now_iso(),
        )

    def save_hierarchy(self, *, with_lap=True):
        day = self.services.save_day(self.make_day())
        setup = self.services.save_setup(self.make_setup(day.id))
        event = self.services.save_event_layout(self.make_event())
        lap = self.services.save_lap(self.make_lap(setup.id, event.id)) if with_lap else None
        return day, setup, event, lap

    @staticmethod
    def set_field(widget, value):
        if isinstance(widget, tk.Text):
            widget.delete("1.0", "end")
            widget.insert("1.0", value)
        elif isinstance(widget, ttk.Combobox):
            widget.set(value)
        else:
            widget.delete(0, "end")
            widget.insert(0, value)

    def select_tree_item(self, tree, item_id):
        tree.selection_set(item_id)
        self.root.update()

    def test_add_day_setup_lap_callbacks_save_the_three_level_hierarchy(self):
        event = self.services.save_event_layout(self.make_event())
        window = self.require_window()

        for action in ("add_day", "add_setup", "add_lap", "save"):
            self.assertIn(action, window.actions)

        window.actions["add_day"].invoke()
        self.set_field(window.fields["date"], "2026-10-05")
        self.set_field(window.fields["location"], "Test Site")
        window.actions["save"].invoke()
        day = self.services.list_days()[0]

        self.select_tree_item(window.tree, f"day:{day.id}")
        window.actions["add_setup"].invoke()
        self.set_field(window.fields["name"], "Baseline")
        self.set_field(window.fields["setup_code"], "A01")
        self.set_field(window.fields["settings_text"], "Cold pressures: 20 psi")
        window.actions["save"].invoke()
        setup = self.services.list_setups(day.id)[0]

        self.select_tree_item(window.tree, f"setup:{setup.id}")
        window.actions["add_lap"].invoke()
        self.assertEqual(window.fields["lap_time"].get(), "")
        self.assertEqual(window.fields["status"].get(), "")
        self.assertTrue(window.actions["move_lap_up"].instate(["disabled"]))
        self.assertTrue(window.actions["move_lap_down"].instate(["disabled"]))
        self.set_field(window.fields["lap_time"], "42.318")
        self.set_field(window.fields["event_layout_id"], window.event_choice_label(event.id))
        self.set_field(window.fields["status"], "valid")
        window.actions["save"].invoke()
        lap = self.services.list_laps(setup.id)[0]

        self.assertEqual((lap.setup_id, lap.event_layout_id, lap.time_ms), (setup.id, event.id, 42318))
        self.assertEqual(window.tree.parent(f"setup:{setup.id}"), f"day:{day.id}")
        self.assertEqual(window.tree.parent(f"lap:{lap.id}"), f"setup:{setup.id}")
        self.assertIn("BEST", window.tree.item(f"lap:{lap.id}", "text"))

    def test_canceling_dirty_selection_restores_selection_and_form_values(self):
        day = self.services.save_day(self.make_day())
        setup = self.services.save_setup(self.make_setup(day.id))
        window = self.require_window()
        self.select_tree_item(window.tree, f"day:{day.id}")
        self.set_field(window.fields["location"], "Changed but unsaved")

        with patch("tkinter.messagebox.askyesnocancel", return_value=None):
            window.tree.selection_set(f"setup:{setup.id}")
            self.root.update()

        self.assertEqual(window.tree.selection(), (f"day:{day.id}",))
        self.assertEqual(window.fields["location"].get(), "Changed but unsaved")
        self.assertEqual(self.services.get_day(day.id).location, "Test Site")

    def test_editing_historical_lap_retains_its_archived_event(self):
        day = self.services.save_day(self.make_day())
        setup = self.services.save_setup(self.make_setup(day.id))
        archived_event = self.services.save_event_layout(self.make_event(layout="Old Layout"))
        other_archived = self.services.save_event_layout(
            self.make_event(layout="Other Old Layout", event_name="Earlier Day")
        )
        lap = self.services.save_lap(self.make_lap(setup.id, archived_event.id))
        self.services.archive_event_layout(archived_event.id)
        self.services.archive_event_layout(other_archived.id)
        window = self.require_window()
        self.select_tree_item(window.tree, f"lap:{lap.id}")

        selected_label = window.event_choice_label(archived_event.id)
        other_label = window.event_choice_label(other_archived.id)
        self.assertTrue(selected_label)
        self.assertNotIn(other_label, window.fields["event_layout_id"]["values"])
        self.assertEqual(window.fields["event_layout_id"].get(), selected_label)

        window.fields["status"].set("invalid")
        window.actions["save"].invoke()

        saved = self.services.get_lap(lap.id)
        self.assertEqual(saved.status, "invalid")
        self.assertEqual(saved.event_layout_id, archived_event.id)
        self.assertIn("INVALID", window.tree.item(f"lap:{lap.id}", "text"))
        self.assertNotIn("BEST", window.tree.item(f"lap:{lap.id}", "text"))

    def test_failed_record_transaction_retains_fields_and_requires_reattachment(self):
        day = self.services.save_day(self.make_day())
        setup = self.services.save_setup(self.make_setup(day.id))
        stale_event = self.services.save_event_layout(self.make_event(layout="Removed Before Save"))
        window = self.require_window()
        self.select_tree_item(window.tree, f"setup:{setup.id}")
        window.actions["add_lap"].invoke()
        self.assertEqual(window.fields["lap_time"].get(), "")
        self.assertEqual(window.fields["status"].get(), "")
        self.assertTrue(window.actions["move_lap_up"].instate(["disabled"]))
        self.assertTrue(window.actions["move_lap_down"].instate(["disabled"]))
        self.set_field(window.fields["lap_time"], "42.318")
        self.set_field(window.fields["event_layout_id"], window.event_choice_label(stale_event.id))
        self.set_field(window.fields["status"], "valid")

        source = self.root_path / "lap-notes.txt"
        source.write_text("source remains", encoding="utf-8")
        with patch("tkinter.filedialog.askopenfilename", return_value=str(source)):
            window.actions["attach_file"].invoke()
        staged = window.staged_attachments[0]
        stored_copy = self.paths.root / staged.relative_path
        self.assertTrue(stored_copy.is_file())

        self.services.delete_event_layout(stale_event.id)
        window.actions["save"].invoke()
        self.assertIn("Select an existing event/layout.", window._field_errors["event_layout_id"].cget("text"))
        self.assertTrue(window.fields["lap_time"].get())
        self.assertIn("reattach", window.status_label.cget("text").lower())
        self.assertEqual(window.staged_attachments, [])
        self.assertFalse(stored_copy.exists())
        self.assertEqual(self.services.list_laps(setup.id), [])

        replacement_event = self.services.save_event_layout(self.make_event(layout="Replacement"))
        window.refresh_tree()
        self.assertEqual(window.fields["lap_time"].get(), "42.318")
        self.set_field(window.fields["event_layout_id"], window.event_choice_label(replacement_event.id))
        with patch("tkinter.filedialog.askopenfilename", return_value=str(source)):
            window.actions["attach_file"].invoke()
        window.actions["save"].invoke()

        saved_lap = self.services.list_laps(setup.id)[0]
        attachments = self.services.list_attachments("lap", saved_lap.id)
        self.assertEqual(saved_lap.event_layout_id, replacement_event.id)
        self.assertEqual([attachment.original_name for attachment in attachments], ["lap-notes.txt"])
        self.assertTrue(self.services.resolve_attachment(attachments[0]).is_file())
        self.assertEqual(source.read_text(encoding="utf-8"), "source remains")

    def test_event_library_creates_duplicates_and_freezes_referenced_history(self):
        day, setup, event, lap = self.save_hierarchy()
        window = self.require_window()
        self.select_tree_item(window.tree, f"lap:{lap.id}")
        self.assertIn("manage_events", window.actions)
        window.actions["manage_events"].invoke()
        manager = window.event_manager
        self.managers.append(manager)
        self.root.update()

        event_item = f"event:{event.id}"
        self.assertTrue(manager.tree.exists(event_item))
        self.select_tree_item(manager.tree, event_item)
        self.assertIn("disabled", manager.fields["track_name"].state())
        self.assertIn("disabled", manager.fields["layout_name"].state())
        self.assertTrue(manager.actions["save"].instate(["disabled"]))
        self.assertFalse(manager.actions["duplicate"].instate(["disabled"]))
        self.assertFalse(manager.actions["archive"].instate(["disabled"]))

        manager.actions["duplicate"].invoke()
        duplicates = [item for item in self.services.list_event_layouts() if item.id != event.id]
        self.assertEqual(len(duplicates), 1)
        duplicate = duplicates[0]
        self.assertFalse(duplicate.archived)
        self.assertNotEqual(duplicate.id, event.id)

        self.select_tree_item(manager.tree, event_item)
        manager.actions["archive"].invoke()
        self.assertTrue(self.services.get_event_layout(event.id).archived)
        self.assertEqual(self.services.get_lap(lap.id).event_layout_id, event.id)
        self.assertIn("Short Loop", window.tree.item(f"lap:{lap.id}", "text"))
        self.assertTrue(manager.tree.exists(event_item))

    def test_attach_and_open_lap_file_uses_a_copied_constrained_attachment(self):
        day, setup, event, lap = self.save_hierarchy()
        window = self.require_window()
        self.select_tree_item(window.tree, f"lap:{lap.id}")

        source = self.root_path / "setup-note.txt"
        source.write_text("driver notes", encoding="utf-8")
        with patch("tkinter.filedialog.askopenfilename", return_value=str(source)):
            window.actions["attach_file"].invoke()
        self.assertEqual(len(window.staged_attachments), 1)
        window.actions["save"].invoke()

        attachment = self.services.list_attachments("lap", lap.id)[0]
        copied_path = self.services.resolve_attachment(attachment)
        self.assertEqual(copied_path.read_text(encoding="utf-8"), "driver notes")
        self.assertEqual(source.read_text(encoding="utf-8"), "driver notes")
        self.assertGreater(window.attachment_list.size(), 0)
        window.attachment_list.selection_set(0)
        with patch("vd_test_log.ui.os.startfile") as startfile:
            window.actions["open_attachment"].invoke()
        startfile.assert_called_once_with(str(copied_path))

    def test_delete_day_confirmation_names_cascade_and_removes_records(self):
        day = self.services.save_day(self.make_day())
        source = self.root_path / "day-note.txt"
        source.write_text("day notes", encoding="utf-8")
        staged = self.services.stage_attachment(source, "file")
        self.services.save_day(day, (staged,))
        setup = self.services.save_setup(self.make_setup(day.id))
        event = self.services.save_event_layout(self.make_event())
        lap = self.services.save_lap(self.make_lap(setup.id, event.id))
        day_attachment = self.services.list_attachments("day", day.id)[0]
        copied_path = self.services.resolve_attachment(day_attachment)
        window = self.require_window()
        self.select_tree_item(window.tree, f"day:{day.id}")

        with patch("tkinter.messagebox.askyesno", return_value=True) as confirm:
            window.actions["delete"].invoke()

        self.assertTrue(confirm.called)
        prompt = " ".join(str(part) for part in confirm.call_args.args).lower()
        self.assertIn("test day", prompt)
        self.assertIn("setup", prompt)
        self.assertIn("lap", prompt)
        self.assertIsNone(self.services.get_day(day.id))
        self.assertIsNone(self.services.get_setup(setup.id))
        self.assertIsNone(self.services.get_lap(lap.id))
        self.assertFalse(copied_path.exists())

    def test_backup_and_csv_buttons_write_files_from_selected_setup(self):
        day, setup, event, lap = self.save_hierarchy()
        window = self.require_window()
        self.select_tree_item(window.tree, f"setup:{setup.id}")

        destination = self.root_path / "test laps.csv"
        with patch("tkinter.filedialog.asksaveasfilename", return_value=str(destination)):
            window.actions["export_csv"].invoke()
        self.assertTrue(destination.is_file())
        self.assertIn("Short Loop", destination.read_text(encoding="utf-8-sig"))

        window.actions["backup"].invoke()
        backups = list(self.paths.backup_root.iterdir())
        self.assertEqual(len(backups), 1)
        self.assertTrue((backups[0] / "test_log.sqlite3").is_file())

    def test_lap_move_buttons_reorder_saved_lap_sequences(self):
        day = self.services.save_day(self.make_day())
        setup = self.services.save_setup(self.make_setup(day.id))
        event = self.services.save_event_layout(self.make_event())
        first = self.services.save_lap(self.make_lap(setup.id, event.id, sequence=1, time_ms=42318))
        second = self.services.save_lap(self.make_lap(setup.id, event.id, sequence=2, time_ms=43111))
        window = self.require_window()
        self.select_tree_item(window.tree, f"lap:{second.id}")

        window.actions["move_lap_up"].invoke()

        laps = self.services.list_laps(setup.id)
        self.assertEqual([lap.id for lap in laps], [second.id, first.id])
        self.assertEqual([lap.sequence for lap in laps], [1, 2])
        self.assertEqual(window.tree.parent(f"lap:{second.id}"), f"setup:{setup.id}")

    def test_event_map_stages_and_opens_through_library_callbacks(self):
        window = self.require_window()
        window.actions["manage_events"].invoke()
        manager = window.event_manager
        self.managers.append(manager)

        manager.actions["new_event"].invoke()
        self.set_field(manager.fields["track_name"], "Test Site")
        self.set_field(manager.fields["layout_name"], "North Loop")
        self.set_field(manager.fields["event_name"], "October Test")
        manager.fields["length_m"].insert(0, "1500")
        source = self.root_path / "north-loop.map"
        source.write_text("course map", encoding="utf-8")
        with patch("tkinter.filedialog.askopenfilename", return_value=str(source)):
            manager.actions["attach_map"].invoke()
        manager.actions["save"].invoke()
        self.set_field(manager.fields["layout_name"], "North Loop Updated")
        manager.actions["save"].invoke()

        event = self.services.list_event_layouts()[0]
        self.assertEqual(event.layout_name, "North Loop Updated")
        attachment = self.services.list_attachments("event_layout", event.id)[0]
        self.assertEqual(attachment.role, "map")
        self.assertEqual(self.services.resolve_attachment(attachment).read_text(encoding="utf-8"), "course map")
        manager.attachment_list.selection_set(0)
        with patch("vd_test_log.event_ui.os.startfile") as startfile:
            manager.actions["open_map"].invoke()
        startfile.assert_called_once_with(str(self.services.resolve_attachment(attachment)))



    def test_duplicate_setup_button_uses_the_service_duplicate(self):
        day = self.services.save_day(self.make_day())
        setup = self.services.save_setup(self.make_setup(day.id))
        window = self.require_window()
        self.select_tree_item(window.tree, f"setup:{setup.id}")

        window.actions["duplicate_setup"].invoke()

        duplicates = self.services.list_setups(day.id)
        self.assertEqual(len(duplicates), 2)
        duplicate = next(item for item in duplicates if item.id != setup.id)
        self.assertEqual(duplicate.name, "Baseline (copy)")
        self.assertIsNone(duplicate.setup_code)
        self.assertEqual(duplicate.settings_text, setup.settings_text)
        self.assertEqual(duplicate.order, 2)
        self.assertTrue(window.tree.exists(f"setup:{duplicate.id}"))

    def test_identical_event_labels_keep_distinct_internal_ids(self):
        day = self.services.save_day(self.make_day())
        setup = self.services.save_setup(self.make_setup(day.id))
        first = self.make_event()
        first = first.__class__(
            id="a1b2c30123456789abcdef0123456789",
            track_name=first.track_name,
            layout_name=first.layout_name,
            event_name=first.event_name,
            event_type=first.event_type,
            length_m=first.length_m,
            notes=first.notes,
            archived=False,
            created_at=first.created_at,
        )
        second = self.make_event()
        second = second.__class__(
            id="a1b2c31123456789abcdef0123456789",
            track_name=second.track_name,
            layout_name=second.layout_name,
            event_name=second.event_name,
            event_type=second.event_type,
            length_m=second.length_m,
            notes=second.notes,
            archived=False,
            created_at=second.created_at,
        )
        self.services.save_event_layout(first)
        self.services.save_event_layout(second)
        window = self.require_window()
        self.select_tree_item(window.tree, f"setup:{setup.id}")
        window.actions["add_lap"].invoke()

        values = tuple(window.fields["event_layout_id"]["values"])
        first_label = window.event_choice_label(first.id)
        second_label = window.event_choice_label(second.id)
        self.assertEqual(len(values), 2)
        self.assertEqual(len(set(values)), 2)
        self.assertNotEqual(first_label, second_label)
        self.set_field(window.fields["event_layout_id"], first_label)
        self.set_field(window.fields["lap_time"], "42.318")
        self.set_field(window.fields["status"], "valid")
        window.actions["save"].invoke()
        self.assertEqual(self.services.list_laps(setup.id)[0].event_layout_id, first.id)

    def test_canceling_close_keeps_dirty_editor_open(self):
        day = self.services.save_day(self.make_day())
        window = self.require_window()
        self.select_tree_item(window.tree, f"day:{day.id}")
        self.set_field(window.fields["location"], "Changed but unsaved")

        with patch("tkinter.messagebox.askyesnocancel", return_value=None):
            self.assertFalse(window.close_request())

        self.assertTrue(self.root.winfo_exists())
        self.assertEqual(window.fields["location"].get(), "Changed but unsaved")
        self.assertEqual(self.services.get_day(day.id).location, "Test Site")

    def test_event_manager_error_labels_do_not_overlap_next_field(self):
        window = self.require_window()
        window.actions["manage_events"].invoke()
        manager = window.event_manager
        self.managers.append(manager)

        error_cell = manager._errors["track_name"].grid_info()
        next_field_cell = manager.fields["layout_name"].grid_info()
        self.assertNotEqual(
            (error_cell["row"], error_cell["column"]),
            (next_field_cell["row"], next_field_cell["column"]),
        )

    def test_main_window_close_guards_dirty_event_manager(self):
        window = self.require_window()
        window.actions["manage_events"].invoke()
        manager = window.event_manager
        self.managers.append(manager)
        manager.actions["new_event"].invoke()
        self.set_field(manager.fields["track_name"], "Unsaved track")

        with patch("tkinter.messagebox.askyesnocancel", return_value=None):
            self.assertFalse(window.close_request())

        self.assertTrue(self.root.winfo_exists())
        self.assertTrue(manager.winfo_exists())
        self.assertEqual(manager.fields["track_name"].get(), "Unsaved track")


if __name__ == "__main__":
    unittest.main()
