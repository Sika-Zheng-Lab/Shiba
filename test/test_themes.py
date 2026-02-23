"""
Unit tests for src/themes.py
Tests cover: COLOR_SCHEMES, get_active_theme, get_css_variables, list_available_themes
"""

import unittest
import os
import sys

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
import themes


class TestColorSchemes(unittest.TestCase):
    def test_color_schemes_is_dict(self):
        self.assertIsInstance(themes.COLOR_SCHEMES, dict)

    def test_expected_themes_present(self):
        expected = {"modern_blue", "scientific_green", "dark_mode", "academic_purple", "minimal_gray"}
        self.assertEqual(set(themes.COLOR_SCHEMES.keys()), expected)

    def test_each_theme_has_required_keys(self):
        for name, theme in themes.COLOR_SCHEMES.items():
            self.assertIn("name", theme, f"{name} missing 'name'")
            self.assertIn("description", theme, f"{name} missing 'description'")
            self.assertIn("colors", theme, f"{name} missing 'colors'")

    def test_each_theme_has_required_colors(self):
        required_colors = [
            "primary", "primary_dark", "secondary", "accent",
            "success", "warning", "error",
            "background", "surface", "text", "text_light"
        ]
        for name, theme in themes.COLOR_SCHEMES.items():
            for color in required_colors:
                self.assertIn(color, theme["colors"], f"{name} missing color '{color}'")

    def test_color_values_are_hex(self):
        for name, theme in themes.COLOR_SCHEMES.items():
            for color_name, color_value in theme["colors"].items():
                self.assertTrue(
                    color_value.startswith("#"),
                    f"{name}.{color_name} = '{color_value}' is not a hex color"
                )


class TestGetActiveTheme(unittest.TestCase):
    def test_returns_dict(self):
        result = themes.get_active_theme()
        self.assertIsInstance(result, dict)

    def test_returns_default_theme(self):
        result = themes.get_active_theme()
        self.assertEqual(result["name"], "Modern Blue")

    def test_active_theme_in_color_schemes(self):
        self.assertIn(themes.ACTIVE_THEME, themes.COLOR_SCHEMES)


class TestGetCssVariables(unittest.TestCase):
    def test_returns_css_string(self):
        result = themes.get_css_variables()
        self.assertIsInstance(result, str)
        self.assertIn(":root", result)

    def test_contains_css_variables(self):
        result = themes.get_css_variables()
        self.assertIn("--primary:", result)
        self.assertIn("--background:", result)

    def test_specific_theme(self):
        result = themes.get_css_variables("dark_mode")
        self.assertIn("Dark Mode", result)
        self.assertIn("#111827", result)  # dark_mode background

    def test_fallback_to_active_theme(self):
        result = themes.get_css_variables("nonexistent_theme")
        # Should fall back to active theme
        self.assertIn("Modern Blue", result)


class TestListAvailableThemes(unittest.TestCase):
    def test_list_available_themes_prints_output(self):
        """list_available_themes prints to stdout without error."""
        import io
        from contextlib import redirect_stdout
        f = io.StringIO()
        with redirect_stdout(f):
            themes.list_available_themes()
        output = f.getvalue()
        self.assertIn("Available", output)
        self.assertIn("modern_blue", output)


if __name__ == "__main__":
    unittest.main()
