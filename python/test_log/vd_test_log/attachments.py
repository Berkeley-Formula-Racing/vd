"""Safe staging, resolving, and cleanup for app-owned attachment copies."""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path
from typing import Iterable

from .models import Attachment, AttachmentRole, DataPaths, StagedAttachment, new_id


class AttachmentPathError(ValueError):
    """A stored attachment path does not stay inside the managed attachment tree."""


class AttachmentMissingError(FileNotFoundError):
    """A stored attachment is missing or is no longer a regular file."""


class AttachmentError(OSError):
    """An attachment copy or cleanup operation failed."""


class AttachmentStore:
    def __init__(self, paths: DataPaths):
        self.paths = paths

    @property
    def root(self) -> Path:
        return Path(self.paths.root).resolve(strict=False)

    @property
    def attachment_root(self) -> Path:
        root = self.root
        configured = Path(self.paths.attachments)
        if _is_link_or_junction(configured):
            raise AttachmentPathError("The attachments folder cannot be a link or junction.")
        attachment_root = configured.resolve(strict=False)
        try:
            attachment_root.relative_to(root)
        except ValueError as error:
            raise AttachmentPathError(
                "The attachments folder must remain inside the selected data folder."
            ) from error
        return attachment_root

    def stage(self, source: Path, role: AttachmentRole) -> StagedAttachment:
        """Copy a selected source into generated app-owned storage."""
        if role not in ("file", "map"):
            raise ValueError("Attachment role must be 'file' or 'map'.")
        source_path = Path(source).expanduser()
        original_name = source_path.name
        if not original_name:
            raise AttachmentError("Choose a file with a name to attach.")
        extension = source_path.suffix
        if not re.fullmatch(r"\.[A-Za-z0-9]{1,16}", extension):
            extension = ""
        attachment_root = self.attachment_root
        attachment_root.mkdir(parents=True, exist_ok=True)
        relative_directory = attachment_root.relative_to(self.root).as_posix()
        relative_path = f"{relative_directory}/{new_id()}{extension}"
        destination = self.resolve_relative(relative_path, require_exists=False)
        created_destination = False
        try:
            with source_path.open("rb") as input_stream:
                with destination.open("xb") as output_stream:
                    created_destination = True
                    shutil.copyfileobj(input_stream, output_stream)
        except OSError:
            if created_destination:
                try:
                    destination.unlink(missing_ok=True)
                except OSError:
                    pass
            raise
        return StagedAttachment(
            role=role,
            original_name=original_name,
            relative_path=relative_path,
        )

    def resolve(self, attachment: Attachment) -> Path:
        """Return an existing file only after checking its resolved path."""
        return self.resolve_relative(attachment.relative_path, require_exists=True)

    def resolve_staged(self, attachment: StagedAttachment) -> Path:
        return self.resolve_relative(attachment.relative_path, require_exists=True)

    def resolve_relative(self, relative_path: str, *, require_exists: bool) -> Path:
        if not isinstance(relative_path, str) or not relative_path.strip():
            raise AttachmentPathError("The stored attachment path is empty.")
        if "\0" in relative_path:
            raise AttachmentPathError("The stored attachment path is invalid.")
        relative = Path(relative_path)
        if relative.is_absolute() or relative.drive:
            raise AttachmentPathError("Stored attachment paths must be relative to the data folder.")
        raw_parts = re.split(r"[\\/]", relative_path)
        if any(part in ("", ".", "..") or ":" in part for part in raw_parts):
            raise AttachmentPathError("Stored attachment paths cannot contain traversal segments.")
        root = self.root
        base = self.attachment_root
        candidate = root
        for part in relative.parts:
            candidate = candidate / part
            if _is_link_or_junction(candidate):
                raise AttachmentPathError("Stored attachment paths cannot pass through links or junctions.")
        try:
            candidate = candidate.resolve(strict=False)
            candidate.relative_to(base)
        except (OSError, ValueError) as error:
            raise AttachmentPathError(
                "The stored attachment path escapes the managed attachments folder."
            ) from error
        if candidate == base:
            raise AttachmentPathError("The stored attachment path must name a file.")
        if require_exists and (not candidate.exists() or not candidate.is_file()):
            raise AttachmentMissingError(f"The attached file is missing: {relative_path}")
        return candidate

    def remove_staged(self, attachments: Iterable[StagedAttachment]) -> list[Path]:
        """Delete staged copies only; the caller must first verify ownership."""
        failures: list[Path] = []
        for attachment in attachments:
            try:
                path = self.resolve_relative(attachment.relative_path, require_exists=False)
            except AttachmentPathError:
                failures.append(Path(attachment.relative_path))
                continue
            try:
                path.unlink(missing_ok=True)
            except OSError:
                failures.append(path)
        return failures

    def remove_committed(self, attachments: Iterable[Attachment]) -> list[Path]:
        """Delete database-referenced copies after their rows have committed away."""
        failures: list[Path] = []
        for attachment in attachments:
            try:
                path = self.resolve_relative(attachment.relative_path, require_exists=False)
            except AttachmentPathError:
                failures.append(Path(attachment.relative_path))
                continue
            try:
                path.unlink(missing_ok=True)
            except OSError:
                failures.append(path)
        return failures


def _is_link_or_junction(path: Path) -> bool:
    if path.is_symlink():
        return True
    is_junction = getattr(os.path, "isjunction", None)
    return bool(is_junction and is_junction(path))
