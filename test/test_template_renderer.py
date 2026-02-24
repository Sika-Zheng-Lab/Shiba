"""
Unit tests for src/template_renderer.py
Tests cover: HTMLTemplateRenderer methods, get_splicing_event_config
"""

import unittest
import os
import sys
import tempfile
import shutil

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
from template_renderer import HTMLTemplateRenderer, get_splicing_event_config


class TestHTMLTemplateRenderer(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a minimal template
        with open(os.path.join(self.tmpdir, "test.html"), "w") as f:
            f.write("<div>{section_id}</div>")
        with open(os.path.join(self.tmpdir, "splicing_section.html"), "w") as f:
            f.write('<section id="{section_id}">{section_title} {section_description} {icon_class} {volcano_content} {scatter_content}</section>')
        with open(os.path.join(self.tmpdir, "summary.html"), "w") as f:
            f.write('<html>{shiba_command}{pca_tpm_content}{pca_psi_content}{splicing_sections}</html>')
        self.renderer = HTMLTemplateRenderer(self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_load_template(self):
        result = self.renderer.load_template("test.html")
        self.assertEqual(result, "<div>{section_id}</div>")

    def test_load_template_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            self.renderer.load_template("nonexistent.html")

    def test_escape_and_format_iframe_content(self):
        lines = ['<div class="test">', '  <p>Hello & World</p>', '</div>']
        result = self.renderer.escape_and_format_iframe_content(lines)
        self.assertIn("&lt;div", result)
        self.assertIn("&amp;", result)

    def test_render_splicing_section(self):
        event_data = {
            "id": "se",
            "icon": "fas fa-dna",
            "title": "Skipped Exon",
            "description": "SE events",
            "volcano_content": "<div>volcano</div>",
            "scatter_content": "<div>scatter</div>",
        }
        result = self.renderer.render_splicing_section(event_data)
        self.assertIn('id="se"', result)
        self.assertIn("Skipped Exon", result)
        self.assertIn("SE events", result)

    def test_render_summary_html(self):
        data = {
            "shiba_command": "shiba run",
            "pca_tpm_content": "<div>pca_tpm</div>",
            "pca_psi_content": "<div>pca_psi</div>",
            "splicing_events": [
                {
                    "id": "se",
                    "icon": "fas fa-dna",
                    "title": "SE",
                    "description": "desc",
                    "volcano_content": "vol",
                    "scatter_content": "scat",
                }
            ],
        }
        result = self.renderer.render_summary_html(data)
        self.assertIn("shiba run", result)
        self.assertIn("pca_tpm", result)

    def test_render_splicing_section_standalone(self):
        event_data = {
            "id": "five",
            "icon": "fas fa-cut",
            "title": "Five Prime",
            "description": "5' splice site",
            "volcano_content": "<div>volcano</div>",
            "scatter_content": "<div>scatter</div>",
        }
        result = self.renderer.render_splicing_section_standalone(event_data)
        self.assertIn('id="five"', result)
        self.assertIn("Five Prime", result)
        self.assertIn("Volcano Plot", result)


class TestGetSplicingEventConfig(unittest.TestCase):
    def test_returns_list(self):
        result = get_splicing_event_config()
        self.assertIsInstance(result, list)

    def test_has_expected_event_types(self):
        result = get_splicing_event_config()
        ids = [e["id"] for e in result]
        for expected in ["se", "five", "three", "mxe"]:
            self.assertIn(expected, ids)

    def test_each_event_has_required_keys(self):
        result = get_splicing_event_config()
        required_keys = {"id", "title", "description", "icon", "code"}
        for event in result:
            for key in required_keys:
                self.assertIn(key, event, f"Event {event.get('id', '?')} missing key '{key}'")


if __name__ == "__main__":
    unittest.main()
