import unittest

from vd_test_log.time_value import format_lap_time, parse_lap_time
from vd_test_log.validation import ValidationError


class LapTimeTests(unittest.TestCase):
    def test_rounds_arbitrary_decimal_seconds_half_up(self):
        self.assertEqual(parse_lap_time("42.3185"), 42319)
        self.assertEqual(
            parse_lap_time("42.3184999999999999999999999999999999999"),
            42318,
        )
        self.assertEqual(parse_lap_time("0.0005"), 1)

    def test_parses_colon_notation_and_canonical_format(self):
        self.assertEqual(parse_lap_time("1:42.318"), 102318)
        self.assertEqual(parse_lap_time("0:59.999"), 59999)
        self.assertEqual(format_lap_time(42318), "0:42.318")
        self.assertEqual(format_lap_time(102318), "1:42.318")
        self.assertEqual(format_lap_time(3_600_000), "60:00.000")

    def test_rejects_zero_rounding_negative_blank_and_nonfinite_values(self):
        for text in ("0", "0.0004", "-1", "", "NaN", "Infinity", "-Infinity"):
            with self.subTest(text=text):
                with self.assertRaises(ValidationError) as raised:
                    parse_lap_time(text)
                self.assertEqual(raised.exception.field, "lap_time")
                self.assertTrue(raised.exception.message)

    def test_rejects_malformed_input_and_colon_seconds_at_or_above_sixty(self):
        for text in ("abc", "1:60.000", "1:99", "1:2.0000", "1:"):
            with self.subTest(text=text):
                with self.assertRaises(ValidationError) as raised:
                    parse_lap_time(text)
                self.assertEqual(raised.exception.field, "lap_time")

    def test_formatter_rejects_nonpositive_milliseconds(self):
        for value in (0, -1):
            with self.subTest(value=value):
                with self.assertRaises(ValidationError) as raised:
                    format_lap_time(value)
                self.assertEqual(raised.exception.field, "lap_time")


if __name__ == "__main__":
    unittest.main()
