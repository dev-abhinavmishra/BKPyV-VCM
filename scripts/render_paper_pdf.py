"""Render docs/CLINICAL_PAPER.md to a typeset PDF (outputs/CLINICAL_PAPER.pdf).

Parses the markdown source into a properly formatted clinical paper:
title block, structured abstract, IMRaD sections, numbered table,
Vancouver-style references, embedded figures, running header/footer.

Requires fpdf2 (pip install fpdf2) plus Liberation Serif TTFs (standard
on Ubuntu images). Re-run after editing the paper:

    .venv/bin/python scripts/render_paper_pdf.py
"""

from __future__ import annotations

import re
from pathlib import Path

from fpdf import FPDF

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "docs" / "CLINICAL_PAPER.md"
OUT = REPO / "outputs" / "CLINICAL_PAPER.pdf"
FIGDIR = REPO / "outputs" / "figures"

FONTS = {
    "": "/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf",
    "B": "/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf",
    "I": "/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf",
    "BI": "/usr/share/fonts/truetype/liberation/LiberationSerif-BoldItalic.ttf",
}
MONO = {
    "": "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf",
    "B": "/usr/share/fonts/truetype/liberation/LiberationMono-Bold.ttf",
}

# Glyphs outside Liberation Serif coverage → plain-text equivalents.
SUBS = "₀₁₂₃₄₅₆₇₈₉"
SUPS = "⁰¹²³⁴⁵⁶⁷⁸⁹⁺"
GLYPH_MAP = {**{c: str(i) for i, c in enumerate(SUBS)},
             **{c: "^" + str(i) for i, c in enumerate(SUPS[:10])},
             "⁺": "^+", "−": "-", "→": "→", "≈": "~", "±": "±"}


def sanitize(text: str) -> str:
    for k, v in GLYPH_MAP.items():
        text = text.replace(k, v)
    return text


# fpdf2's markdown=True handles **bold** but leaves *italic* markers as
# literal asterisks — strip single-star emphasis before rendering.
ITALIC_RE = re.compile(r"(?<!\*)\*([^*\n]+)\*(?!\*)")


def md(text: str) -> str:
    return ITALIC_RE.sub(r"\1", text)


class Paper(FPDF):
    def __init__(self) -> None:
        super().__init__(format="letter", unit="mm")
        self.set_auto_page_break(True, margin=20)
        self.set_margins(22, 22, 22)
        for style, path in FONTS.items():
            self.add_font("serif", style, path)
        for style, path in MONO.items():
            self.add_font("mono", style, path)
        self.set_font("serif", "", 10)

    def header(self) -> None:
        if self.page_no() == 1:
            return
        self.set_font("serif", "I", 8)
        self.set_text_color(110, 105, 96)
        self.cell(0, 4, "BKPyV-VCM  ·  Clinical paper  ·  v2.1", align="L")
        self.cell(0, 4, "Mishra — October 2026", align="R", new_x="LMARGIN", new_y="NEXT")
        self.set_draw_color(190, 185, 172)
        self.line(22, 12, 188, 12)
        self.ln(4)

    def footer(self) -> None:
        self.set_y(-14)
        self.set_font("serif", "", 8)
        self.set_text_color(110, 105, 96)
        self.cell(0, 8, f"Page {self.page_no()} of {{nb}}", align="C")
        self.set_text_color(30, 30, 28)

    def h2(self, text: str) -> None:
        self.ln(2)
        self.set_font("serif", "B", 13)
        self.set_text_color(13, 76, 71)
        self.multi_cell(0, 6, text, new_x="LMARGIN", new_y="NEXT")
        self.set_text_color(30, 30, 28)
        self.ln(0.5)

    def h3(self, text: str) -> None:
        self.ln(1.5)
        self.set_font("serif", "B", 11)
        self.multi_cell(0, 5.5, text, new_x="LMARGIN", new_y="NEXT")
        self.ln(0.4)

    def h4(self, text: str) -> None:
        self.set_font("serif", "BI", 10)
        self.multi_cell(0, 5.2, text, new_x="LMARGIN", new_y="NEXT")

    def body(self, text: str, indent: float = 0, size: float = 10,
             style: str = "") -> None:
        self.set_font("serif", style, size)
        if indent:
            self.set_x(self.l_margin + indent)
        self.multi_cell(self.w - self.l_margin - self.r_margin - indent,
                        4.6, md(text), markdown=True, align="J", new_x="LMARGIN", new_y="NEXT")
        self.ln(1.2)

    def bullet(self, text: str) -> None:
        x = self.l_margin
        self.set_font("serif", "", 10)
        self.set_x(x)
        self.cell(5, 4.6, "•")
        self.multi_cell(self.w - self.l_margin - self.r_margin - 5,
                        4.6, md(text), markdown=True, align="J", new_x="LMARGIN", new_y="NEXT")
        self.ln(0.8)


def render_table(pdf: Paper, rows: list[list[str]]) -> None:
    n = len(rows[0])
    if n == 6:
        widths = [22, 26, 24, 32, 30, 32]
    elif n == 5:
        widths = [30, 40, 34, 36, 26]
    else:
        widths = [166 / n] * n
    pdf.set_font("serif", "B", 9)
    pdf.set_fill_color(238, 234, 224)
    pdf.set_draw_color(150, 145, 130)
    for r, row in enumerate(rows):
        if r == 0:
            pdf.set_font("serif", "B", 9)
        elif r == len(rows) - 1:
            pdf.set_font("serif", "B", 9)
        else:
            pdf.set_font("serif", "", 9)
        for c, cell in enumerate(row):
            fill = r == 0 or r == len(rows) - 1 or (r % 2 == 0)
            pdf.cell(widths[c], 6.5, re.sub(r"\*\*", "", cell), border=1,
                     align="C", fill=fill)
        pdf.ln()
    pdf.ln(2)


def main() -> None:
    lines = sanitize(MD.read_text()).splitlines()
    pdf = Paper()
    pdf.alias_nb_pages()
    pdf.add_page()

    i = 0
    # Title block: first line is "# <title>"
    title = lines[0].lstrip("# ").strip()
    pdf.set_font("serif", "B", 16.5)
    pdf.multi_cell(0, 7.5, title, align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(1)
    pdf.set_font("serif", "I", 10)
    pdf.multi_cell(0, 5, "Abhinav Mishra", align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("serif", "", 9)
    pdf.set_text_color(90, 86, 78)
    pdf.multi_cell(0, 4.6,
                   "Health Careers High School, San Antonio, Texas  ·  "
                   "Regeneron ISEF 2027 candidate  ·  October 2026", align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.multi_cell(0, 4.6,
                   "Canonical living document — updated as the project evolves",
                   align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(30, 30, 28)
    pdf.set_draw_color(150, 145, 130)
    pdf.line(40, pdf.get_y() + 2, 170, pdf.get_y() + 2)
    pdf.ln(5)

    i = 1
    fig_counter = 0
    table_rows: list[list[str]] = []
    in_table = False

    while i < len(lines):
        line = lines[i].rstrip()
        stripped = line.strip()

        if in_table:
            if stripped.startswith("|"):
                cells = [c.strip() for c in stripped.strip("|").split("|")]
                if not all(re.fullmatch(r":?-{3,}:?", c) for c in cells):
                    table_rows.append(cells)
                i += 1
                continue
            render_table(pdf, table_rows)
            table_rows = []
            in_table = False
            continue

        if stripped.startswith("|") and "|" in stripped[1:]:
            in_table = True
            table_rows = [[c.strip() for c in stripped.strip("|").split("|")]]
            i += 1
            continue

        if not stripped or stripped == "---":
            i += 1
            continue
        if stripped.startswith("*") and stripped.endswith("*") and len(stripped) < 400 \
                and not stripped.startswith("**"):
            pdf.body(stripped.strip("*"), style="I")
            i += 1
            continue

        m = re.match(r"^(#{1,4})\s+(.*)", stripped)
        if m:
            level, text = len(m.group(1)), sanitize(m.group(2))
            if level == 2:
                pdf.h2(text)
            elif level == 3:
                pdf.h3(text)
            else:
                pdf.h4(text)
            i += 1
            continue

        fig_m = re.search(r"`outputs/figures/([\w.]+)`", stripped)
        if fig_m and fig_counter < 4:
            fig_counter += 1
            png = FIGDIR / fig_m.group(1)
            cap = md(stripped[2:] if stripped.startswith("- ") else stripped)
            cap = re.sub(r"\*\*", "", cap)
            cap = re.sub(r"\([^)]+\)", "", cap).strip()
            if png.exists():
                if pdf.get_y() > 190:
                    pdf.add_page()
                w = 158
                pdf.image(str(png), x=(pdf.w - w) / 2, w=w)
                pdf.ln(1)
            pdf.set_font("serif", "", 8.5)
            pdf.set_text_color(80, 76, 68)
            pdf.set_x((pdf.w - 150) / 2)
            pdf.multi_cell(150, 4, cap, align="C",
                           new_x="LMARGIN", new_y="NEXT")
            pdf.set_text_color(30, 30, 28)
            pdf.ln(3)
            i += 1
            continue

        if stripped.startswith("- "):
            pdf.bullet(stripped[2:])
            i += 1
            continue

        if re.match(r"^\[\d+\]", stripped):
            pdf.body(stripped, indent=4, size=8.5)
            i += 1
            continue

        if stripped.startswith("**Table") or stripped.startswith("Table "):
            pdf.set_font("serif", "B", 9.5)
            pdf.multi_cell(0, 5, md(re.sub(r"\*\*", "", stripped)), align="C",
                           new_x="LMARGIN", new_y="NEXT")
            pdf.ln(1)
            i += 1
            continue

        pdf.body(stripped)
        i += 1

    if in_table and table_rows:
        render_table(pdf, table_rows)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pdf.output(str(OUT))
    print(f"wrote {OUT} ({OUT.stat().st_size//1024} KB, {pdf.page_no()} pages)")


if __name__ == "__main__":
    main()
