"""Exact decimal parsing and formatting for lap times."""

from __future__ import annotations

import re
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP, localcontext

from .validation import ValidationError


_SECONDS_PATTERN = re.compile(r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)\Z")
_COLON_SECONDS_PATTERN = re.compile(r"(?:[0-5]?[0-9])(?:\.[0-9]{1,3})?\Z")


def parse_lap_time(text: str) -> int:
    """Parse seconds or m:ss.sss and return positive rounded milliseconds."""
    if not isinstance(text, str):
        raise ValidationError("lap_time", "Enter seconds or m:ss.sss.")
    value = text.strip()
    if not value:
        raise ValidationError("lap_time", "Enter a lap time.")

    with localcontext() as context:
        context.prec = max(28, len(value) + 8)
        if ":" in value:
            if value.count(":") != 1:
                raise ValidationError("lap_time", "Use m:ss.sss with seconds below 60.")
            minute_text, second_text = value.split(":", 1)
            if not re.fullmatch(r"[0-9]+", minute_text) or not _COLON_SECONDS_PATTERN.fullmatch(second_text):
                raise ValidationError("lap_time", "Use m:ss.sss with seconds below 60.")
            try:
                seconds = Decimal(minute_text) * 60 + Decimal(second_text)
            except InvalidOperation as error:
                raise ValidationError("lap_time", "Enter a valid lap time.") from error
        elif _SECONDS_PATTERN.fullmatch(value):
            try:
                seconds = Decimal(value)
            except InvalidOperation as error:
                raise ValidationError("lap_time", "Enter a valid lap time.") from error
        else:
            try:
                numeric_value = Decimal(value)
            except InvalidOperation as error:
                raise ValidationError("lap_time", "Enter seconds or m:ss.sss.") from error
            if not numeric_value.is_finite() or numeric_value <= 0:
                raise ValidationError("lap_time", "Enter a positive finite lap time.")
            raise ValidationError("lap_time", "Enter seconds or m:ss.sss.")

        if not seconds.is_finite() or seconds <= 0:
            raise ValidationError("lap_time", "Enter a positive finite lap time.")
        try:
            milliseconds = int(
                (seconds * 1000).quantize(Decimal("1"), rounding=ROUND_HALF_UP)
            )
        except (InvalidOperation, OverflowError) as error:
            raise ValidationError("lap_time", "Enter a lap time that can be represented in milliseconds.") from error

    if milliseconds <= 0:
        raise ValidationError("lap_time", "Lap time rounds to zero milliseconds.")
    return milliseconds


def format_lap_time(time_ms: int) -> str:
    """Format positive integer milliseconds as m:ss.sss without float math."""
    if isinstance(time_ms, bool) or not isinstance(time_ms, int) or time_ms <= 0:
        raise ValidationError("lap_time", "Lap time must be positive milliseconds.")
    total_seconds, milliseconds = divmod(time_ms, 1000)
    minutes, seconds = divmod(total_seconds, 60)
    return f"{minutes}:{seconds:02d}.{milliseconds:03d}"
