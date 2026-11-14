"""Styled Excel exports, including header-only sheets and atomic replacement."""
import os
from pathlib import Path
import tempfile


def write_workbook(path, sheets, freeze_panes="B2"):
    import pandas as pd
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter

    path = Path(path)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".xlsx", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with pd.ExcelWriter(temporary, engine="openpyxl") as writer:
            for name, frame in sheets:
                frame.to_excel(writer, index=False, sheet_name=name)
                sheet = writer.sheets[name]
                sheet.freeze_panes = freeze_panes
                for index in range(1, len(frame.columns) + 1):
                    sheet.column_dimensions[get_column_letter(index)].width = 20
                for row in sheet:
                    for cell in row:
                        cell.alignment = Alignment(horizontal="left", wrap_text=False)
                        cell.font = Font(name="Calibri", size=11, bold=cell.row == 1)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
