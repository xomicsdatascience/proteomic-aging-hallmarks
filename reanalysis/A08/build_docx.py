import os
import re
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

md_path = "reanalysis/A08/Aging_Proteomics_Manuscript_Revised.md"
docx_path = "reanalysis/A08/Aging_Proteomics_Manuscript_Revised.docx"

with open(md_path, "r", encoding="utf-8") as f:
    text = f.read()

assert "\u2014" not in text, "Em dash found in manuscript markdown!"

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Calibri"
style.font.size = Pt(11)

for line in text.splitlines():
    s = line.strip()
    if not s or s == "---":
        continue
    if s.startswith("# "):
        p = doc.add_heading(s[2:], level=0)
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    elif s.startswith("## "):
        doc.add_heading(s[3:], level=1)
    elif s.startswith("### "):
        doc.add_heading(s[4:], level=2)
    elif re.match(r"^\d+\.\s+", s) or s.startswith("- "):
        clean = re.sub(r"\*\*(.*?)\*\*", r"\1", s)
        clean = clean.replace("`", "")
        if s.startswith("- "):
            doc.add_paragraph(clean[2:], style="List Bullet")
        else:
            doc.add_paragraph(clean)
    else:
        clean = re.sub(r"\*\*(.*?)\*\*", r"\1", s)
        clean = clean.replace("`", "")
        doc.add_paragraph(clean)

doc.add_page_break()
doc.add_heading("Figures", level=1)
figs = [
    ("Figure 1", "reanalysis/A08/figures/fig1_provenance_and_qc.png"),
    ("Figure 2", "reanalysis/A08/figures/fig2_evidence_and_coverage.png"),
    ("Figure 3", "reanalysis/A08/figures/fig3_compound_response_atlas.png"),
    ("Figure 4", "reanalysis/A08/figures/fig4_independent_aging_benchmarks.png"),
    ("Figure 5", "reanalysis/A08/figures/fig5_candidate_mechanisms.png"),
    ("Supplementary Figure S1", "reanalysis/A08/figures/figS1_full_161_compound_atlas.png"),
    ("Supplementary Figure S2", "reanalysis/A08/figures/figS2_sensitivity_and_composition.png"),
]
for title, fpath in figs:
    if os.path.exists(fpath):
        doc.add_heading(title, level=2)
        doc.add_picture(fpath, width=Inches(6.4))

doc.save(docx_path)
print(f"Saved {docx_path} ({os.path.getsize(docx_path)} bytes)")
