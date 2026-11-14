"""Workbook contents and failure behavior with empty splicing results."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lib import shibalib
from lib.excel import write_workbook


@unittest.skipUnless(importlib.util.find_spec("openpyxl"), "openpyxl is required")
class TestExcel(unittest.TestCase):
    def test_empty_and_mixed_splicing_sheets(self):
        from openpyxl import load_workbook
        names = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]
        empty = pd.DataFrame(columns=["event_id", "gene_name", "q"])
        filled = pd.DataFrame([["RI_1", "Gene1", 0.01]], columns=empty.columns)
        for save in (shibalib.save_excel, shibalib.save_excel_sc):
            for mixed in (False, True):
                with self.subTest(exporter=save.__name__, mixed=mixed), tempfile.TemporaryDirectory() as tmp:
                    frames = [empty] * 8
                    if mixed:
                        frames[4] = filled
                    save(tmp, *frames)
                    book = load_workbook(Path(tmp) / "results.xlsx")
                    try:
                        self.assertEqual(book.sheetnames, names)
                        for name in names:
                            sheet = book[name]
                            self.assertEqual([c.value for c in sheet[1]], list(empty.columns))
                            self.assertEqual(sheet.freeze_panes, "B2")
                            self.assertEqual(sheet.max_row, 2 if mixed and name == "RI" else 1)
                        if mixed:
                            self.assertEqual(book["RI"]["A2"].value, "RI_1")
                            self.assertEqual(book["RI"]["C2"].value, 0.01)
                    finally:
                        book.close()

    def test_failure_does_not_publish_partial_workbook(self):
        original = pd.DataFrame.to_excel
        def fail_second(frame, *args, **kwargs):
            if kwargs["sheet_name"] == "FIVE":
                raise OSError("simulated write failure")
            return original(frame, *args, **kwargs)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "results.xlsx"
            path.write_bytes(b"previous complete workbook")
            with patch.object(pd.DataFrame, "to_excel", new=fail_second):
                with self.assertRaisesRegex(OSError, "simulated write failure"):
                    shibalib.save_excel(tmp, *([pd.DataFrame(columns=["event_id", "q"])] * 8))
            self.assertEqual(path.read_bytes(), b"previous complete workbook")
            self.assertEqual(list(Path(tmp).iterdir()), [path])

    def test_empty_expression_sheets(self):
        from openpyxl import load_workbook
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "expression.xlsx"
            frame = pd.DataFrame(columns=["gene_id", "gene_name", "sample"])
            write_workbook(path, [("TPM", frame), ("CPM", frame)], freeze_panes="C2")
            book = load_workbook(path)
            try:
                self.assertEqual(book.sheetnames, ["TPM", "CPM"])
                self.assertEqual(book["TPM"].freeze_panes, "C2")
            finally:
                book.close()
