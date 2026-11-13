import os
import re
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

md_path = "results/Aging_Hallmarks_Proteomics_Report.md"
docx_path = "results/Aging_Hallmarks_Proteomics_Report.docx"

with open(md_path, "r", encoding="utf-8") as f:
    lines = f.read().splitlines()

doc = Document()

# Set clean academic margins & base font
for section in doc.sections:
    section.top_margin = Inches(0.85)
    section.bottom_margin = Inches(0.85)
    section.left_margin = Inches(0.9)
    section.right_margin = Inches(0.9)

style = doc.styles["Normal"]
style.font.name = "Calibri"
style.font.size = Pt(11)
style.paragraph_format.space_after = Pt(6)
style.paragraph_format.line_spacing = 1.15


def add_formatted_runs(p, text, italic_default=False):
    # Clean LaTeX inline symbols for Word readability
    text = text.replace(r"\log_2", "log2").replace(r"\text{FC}", "FC").replace(r"\text{SD}", "SD")
    text = text.replace(r"\times 10^{-15}", "× 10^-15").replace(r"\times 10^{-14}", "× 10^-14")
    text = text.replace(r"\times 10^{-13}", "× 10^-13").replace(r"\times 10^{-11}", "× 10^-11")
    text = text.replace(r"\times 10^{-10}", "× 10^-10").replace(r"\times 10^{-7}", "× 10^-7")
    text = text.replace(r"\mu - 1.8\sigma", "μ - 1.8σ").replace(r"0.3\sigma", "0.3σ")
    text = text.replace(r"\ge", "≥").replace(r"\approx", "≈").replace(r"\alpha", "α").replace(r"\gamma", "γ")
    text = text.replace("$", "")

    # Tokenize **bold**, *italic*, `code`
    pattern = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)")
    pos = 0
    for m in pattern.finditer(text):
        if m.start() > pos:
            r = p.add_run(text[pos:m.start()])
            r.italic = italic_default
        tok = m.group(0)
        if tok.startswith("**") and tok.endswith("**"):
            r = p.add_run(tok[2:-2])
            r.bold = True
            r.italic = italic_default
        elif tok.startswith("*") and tok.endswith("*"):
            r = p.add_run(tok[1:-1])
            r.italic = True
        elif tok.startswith("`") and tok.endswith("`"):
            r = p.add_run(tok[1:-1])
            r.font.name = "Consolas"
            r.font.size = Pt(10)
            r.italic = italic_default
        pos = m.end()
    if pos < len(text):
        r = p.add_run(text[pos:])
        r.italic = italic_default


i = 0
while i < len(lines):
    line = lines[i].strip()
    if not line or line == "---":
        i += 1
        continue

    # Headings
    if line.startswith("# "):
        p = doc.add_heading(level=0)
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
        r = p.add_run(line[2:].strip())
        r.bold = True
        r.font.size = Pt(17)
        r.font.color.rgb = RGBColor(0x1A, 0x36, 0x5D)
        i += 1
        continue
    if line.startswith("## "):
        h = doc.add_heading(line[3:].strip(), level=1)
        for r in h.runs:
            r.font.color.rgb = RGBColor(0x1A, 0x36, 0x5D)
        i += 1
        continue
    if line.startswith("### "):
        h = doc.add_heading(line[4:].strip(), level=2)
        for r in h.runs:
            r.font.color.rgb = RGBColor(0x2B, 0x4C, 0x7E)
        i += 1
        continue

    # Image + optional caption on next line
    img_m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)$", line)
    if img_m:
        img_file = os.path.join("results", img_m.group(2))
        if os.path.exists(img_file):
            p_img = doc.add_paragraph()
            p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_img.add_run().add_picture(img_file, width=Inches(6.3))
        i += 1
        if i < len(lines) and lines[i].strip().startswith("*") and lines[i].strip().endswith("*"):
            cap_text = lines[i].strip().strip("*")
            p_cap = doc.add_paragraph()
            p_cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_cap.paragraph_format.space_after = Pt(10)
            add_formatted_runs(p_cap, cap_text, italic_default=True)
            for r in p_cap.runs:
                r.font.size = Pt(9.5)
                r.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
            i += 1
        continue

    # Bullet or numbered list
    if line.startswith("- "):
        p = doc.add_paragraph(style="List Bullet")
        add_formatted_runs(p, line[2:].strip())
        i += 1
        continue
    num_m = re.match(r"^(\d+)\.\s+(.*)$", line)
    if num_m:
        p = doc.add_paragraph(style="List Number")
        add_formatted_runs(p, num_m.group(2).strip())
        i += 1
        continue

    # Normal paragraph
    p = doc.add_paragraph()
    add_formatted_runs(p, line)
    i += 1

# Append Supplementary Figure S1 on a new page
doc.add_page_break()
h_sup = doc.add_heading("Supplementary Information", level=1)
for r in h_sup.runs:
    r.font.color.rgb = RGBColor(0x1A, 0x36, 0x5D)

s1_path = "results/figS1_aging_hallmarks_nes_heatmap_full.png"
if os.path.exists(s1_path):
    p_s1 = doc.add_paragraph()
    p_s1.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p_s1.add_run().add_picture(s1_path, width=Inches(5.8))
    p_s1_cap = doc.add_paragraph()
    p_s1_cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_formatted_runs(
        p_s1_cap,
        "Supplementary Figure S1. Complete clustered NES heatmap of all 165 drug treatments with at least one significant pathway alteration (FDR q < 0.05; * q < 0.05, ** q < 0.01).",
        italic_default=True,
    )
    for r in p_s1_cap.runs:
        r.font.size = Pt(9.5)
        r.font.color.rgb = RGBColor(0x44, 0x44, 0x44)

doc.save(docx_path)
print("Saved Word document:", docx_path)
