import os
import re
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

md_path = "reanalysis/A08/Aging_Proteomics_Manuscript_Revised.md"
docx_path = "reanalysis/A08/Aging_Proteomics_Manuscript_Revised.docx"

with open(md_path, "r", encoding="utf-8") as f:
    lines = f.read().splitlines()

assert "\u2014" not in "\n".join(lines), "Em dash found in manuscript markdown!"

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Calibri"
style.font.size = Pt(11)

i = 0
while i < len(lines):
    s = lines[i].strip()
    if not s or s == "---":
        i += 1
        continue
    # Markdown table detection
    if s.startswith("|") and i + 1 < len(lines) and lines[i + 1].strip().startswith("|"):
        table_lines = []
        while i < len(lines) and lines[i].strip().startswith("|"):
            row_str = lines[i].strip()
            if not re.match(r"^\|[\s:\-\|]+\|$", row_str):
                cells = [c.strip() for c in row_str.strip("|").split("|")]
                table_lines.append(cells)
            i += 1
        if table_lines:
            tbl = doc.add_table(rows=len(table_lines), cols=len(table_lines[0]))
            tbl.style = "Table Grid"
            for r_idx, row_cells in enumerate(table_lines):
                for c_idx, val in enumerate(row_cells):
                    clean_val = re.sub(r"\*\*(.*?)\*\*", r"\1", val).replace("`", "")
                    cell = tbl.cell(r_idx, c_idx)
                    cell.text = clean_val
                    if r_idx == 0:
                        for run in cell.paragraphs[0].runs:
                            run.bold = True
        continue
    if s.startswith("# "):
        p = doc.add_heading(s[2:], level=0)
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    elif s.startswith("## "):
        doc.add_heading(s[3:], level=1)
    elif s.startswith("### "):
        doc.add_heading(s[4:], level=2)
    elif s.startswith("- "):
        clean = re.sub(r"\*\*(.*?)\*\*", r"\1", s[2:]).replace("`", "")
        doc.add_paragraph(clean, style="List Bullet")
    else:
        clean = re.sub(r"\*\*(.*?)\*\*", r"\1", s).replace("`", "")
        doc.add_paragraph(clean)
    i += 1

doc.add_page_break()
doc.add_heading("Figures", level=1)
figs = [
    ("Figure 1", "reanalysis/A08/figures/fig1_program_coverage.png"),
    ("Figure 2", "reanalysis/A08/figures/fig2_highlighted_heatmap.png"),
    ("Figure 3", "reanalysis/A08/figures/fig3_nes_distributions.png"),
    ("Figure 4", "reanalysis/A08/figures/fig4_kinase_inhibitor_panel.png"),
    ("Supplementary Figure S1", "reanalysis/A08/figures/figS1_qc_and_null_calibration.png"),
    ("Supplementary Figure S2", "reanalysis/A08/figures/figS2_full_160_compound_atlas.png"),
    ("Supplementary Figure S3", "reanalysis/A08/figures/figS3_benchmarks_and_sensitivity.png"),
]
for title, fpath in figs:
    if os.path.exists(fpath):
        doc.add_heading(title, level=2)
        doc.add_picture(fpath, width=Inches(6.4))

doc.save(docx_path)
print(f"Saved {docx_path} ({os.path.getsize(docx_path)} bytes)")
