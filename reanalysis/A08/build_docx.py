import os
import re
import shutil
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

md_path = "reanalysis/A08/Aging_Proteomics_Manuscript_Revised.md"
docx_path = "reanalysis/A08/Aging_Proteomics_Manuscript_Revised.docx"
base_dir = os.path.dirname(md_path)

with open(md_path, "r", encoding="utf-8") as f:
    raw_text = f.read()

if "\u2014" in raw_text:
    raw_text = raw_text.replace("\u2014", " - ")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(raw_text)

lines = raw_text.splitlines()
assert "\u2014" not in "\n".join(lines), "Em dash found in manuscript markdown!"

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Calibri"
style.font.size = Pt(11)

def add_formatted_paragraph(doc, text, style_name=None, italic_default=False):
    p = doc.add_paragraph(style=style_name) if style_name else doc.add_paragraph()
    # Clean backticks
    text_clean = text.replace("`", "")
    # Split by **bold**
    parts = re.split(r"(\*\*.*?\*\*)", text_clean)
    for part in parts:
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            run = p.add_run(part[2:-2])
            run.bold = True
            if italic_default:
                run.italic = True
        else:
            # Handle *italic*
            subparts = re.split(r"(\*[^\*]+?\*)", part)
            for sp in subparts:
                if not sp:
                    continue
                if sp.startswith("*") and sp.endswith("*") and len(sp) > 2:
                    run = p.add_run(sp[1:-1])
                    run.italic = True
                else:
                    run = p.add_run(sp)
                    if italic_default:
                        run.italic = True
    return p

i = 0
while i < len(lines):
    s = lines[i].strip()
    if not s or s == "---":
        i += 1
        continue

    # Inline Markdown image detection: ![alt](path)
    img_match = re.match(r"^!\[(.*?)\]\((.*?)\)$", s)
    if img_match:
        alt_text, rel_img_path = img_match.groups()
        img_full_path = os.path.join(base_dir, rel_img_path)
        if not os.path.exists(img_full_path):
            img_full_path = rel_img_path
        if os.path.exists(img_full_path):
            p_img = doc.add_paragraph()
            p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
            run_img = p_img.add_run()
            run_img.add_picture(img_full_path, width=Inches(6.4))
        i += 1
        # Check if next non-empty line is the Figure Legend
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i < len(lines) and (lines[i].strip().startswith("**Figure") or lines[i].strip().startswith("**Supplementary Figure")):
            leg_line = lines[i].strip()
            p_leg = add_formatted_paragraph(doc, leg_line)
            p_leg.paragraph_format.space_after = Pt(12)
            i += 1
        continue

    # Display math block $$ ... $$
    if s.startswith("$$") and s.endswith("$$") and len(s) > 4:
        p_math = doc.add_paragraph(s.strip("$").strip())
        p_math.alignment = WD_ALIGN_PARAGRAPH.CENTER
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
                    for run in cell.paragraphs[0].runs:
                        run.font.size = Pt(9.5)
                        if r_idx == 0:
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
        add_formatted_paragraph(doc, s[2:], style_name="List Bullet")
    else:
        add_formatted_paragraph(doc, s)
    i += 1

for out_p in [docx_path, "Aging_Hallmarks_Proteomics_Report.docx", "results/Aging_Hallmarks_Proteomics_Report.docx", "reanalysis/A08/Aging_Proteomics_Manuscript_Final.docx"]:
    try:
        doc.save(out_p)
        print(f"Saved {out_p} ({os.path.getsize(out_p)} bytes) with inline figures and immediate legends.")
    except PermissionError:
        print(f"Note: {out_p} is currently open in Word; saved to reanalysis/A08/Aging_Proteomics_Manuscript_Final.docx instead.")
