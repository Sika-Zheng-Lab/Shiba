"""Styled Excel exports, including header-only sheets and atomic replacement."""
import os
from pathlib import Path
import tempfile


def write_workbook(path, sheets, freeze_panes="B2"):
    from styleframe import StyleFrame, Styler, utils

    path = Path(path)
    style = Styler(horizontal_alignment=utils.horizontal_alignments.left, wrap_text=False)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".xlsx", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with StyleFrame.ExcelWriter(temporary) as writer:
            for name, frame in sheets:
                styled = StyleFrame(frame)
                styled.set_column_width(columns=frame.columns, width=20)
                styled.apply_column_style(cols_to_style=frame.columns, styler_obj=style, style_header=True)
                styled.to_excel(writer, index=False, sheet_name=name)
                # openpyxl permits a freeze boundary beyond the last populated
                # row; StyleFrame's argument rejects header-only sheets.
                writer.sheets[name].freeze_panes = freeze_panes
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
