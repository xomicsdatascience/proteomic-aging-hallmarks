"""Export one canonical manuscript and a readable supplementary figure document."""
import re
import sys
import shutil
from io import BytesIO
from zipfile import ZipFile, ZIP_DEFLATED
from pathlib import Path
from PIL import Image
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from analysis_common import ROOT

BASE=ROOT/'reanalysis/A08'


def package_parts(data):
    with ZipFile(BytesIO(data)) as package:
        return {name:package.read(name) for name in package.namelist()}


def save_document(doc,path,compact_images=False):
    """Preserve unchanged packages; compact main images at native pixel size."""
    buffer=BytesIO();doc.save(buffer);data=buffer.getvalue()
    if compact_images:
        output=BytesIO()
        with ZipFile(BytesIO(data)) as source,ZipFile(output,'w',compression=ZIP_DEFLATED,compresslevel=9) as target:
            for entry in source.infolist():
                part=source.read(entry.filename)
                if entry.filename.startswith('word/media/') and entry.filename.endswith('.png'):
                    with Image.open(BytesIO(part)) as image:
                        indexed=image.convert('RGB').quantize(colors=128,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE)
                    image_bytes=BytesIO();indexed.save(image_bytes,format='PNG',optimize=True);part=image_bytes.getvalue()
                target.writestr(entry,part,compress_type=ZIP_DEFLATED,compresslevel=9)
        data=output.getvalue()
    if path.exists() and package_parts(path.read_bytes())==package_parts(data):
        return
    path.write_bytes(data)


def styled_document():
    doc=Document();sec=doc.sections[0]
    sec.page_width=Inches(8.5);sec.page_height=Inches(11)
    sec.top_margin=sec.bottom_margin=sec.left_margin=sec.right_margin=Inches(.75)
    normal=doc.styles['Normal'];normal.font.name='Calibri';normal.font.size=Pt(11)
    normal.paragraph_format.space_after=Pt(7)
    normal.paragraph_format.line_spacing=1.06
    normal.paragraph_format.keep_together=True
    for style in doc.styles:
        for border in style.element.findall('.//'+qn('w:pBdr')):
            border.getparent().remove(border)
    for name,size in [('Title',17),('Heading 1',13),('Heading 2',11.5),('Heading 3',11)]:
        st=doc.styles[name];st.font.name='Calibri';st.font.size=Pt(size);st.font.color.rgb=RGBColor(0,0,0)
        st.font.underline=False
        if name!='Title':st.font.bold=True
    doc.core_properties.title='Proteomic profiling distinguishes drug responses across aging associated cellular programs'
    doc.core_properties.author='Jesse G. Meyer'
    return doc


def formatted(doc,text,style=None):
    p=doc.add_paragraph(style=style)
    # Only alphabetic bold spans are formatting. Literal * and ** significance markers survive.
    parts=re.split(r'(\*\*(?=[A-Za-z])[^*]+\*\*)',text)
    for part in parts:
        if part.startswith('**') and part.endswith('**') and len(part)>4:
            p.add_run(part[2:-2]).bold=True
        else:p.add_run(part)
    return p


def add_image(doc,path,width=6.9):
    p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next=True;p.paragraph_format.space_after=Pt(5)
    p.add_run().add_picture(str(path),width=Inches(width))
    return p


def main():
    text=(BASE/'Aging_Proteomics_Manuscript_Revised.md').read_text()
    if '$' in text or '\\frac' in text or '\u2014' in text:
        raise ValueError('Unrendered math or prohibited punctuation in manuscript')
    lines=text.splitlines();doc=styled_document();captions={};i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        image=re.match(r'^!\[([^\]]*)\]\(([^)]+)\)$',line)
        if image:
            # Keep each image with its legend and let Word fill the preceding page.
            add_image(doc,BASE/image.group(2),width=6.5 if image.group(1)=='Figure 3' else 6.9);i+=1
            while i<len(lines) and not lines[i].strip():i+=1
            if i<len(lines) and lines[i].startswith('**Figure'):
                formatted(doc,lines[i]);i+=1
            continue
        if line.startswith('# '):doc.add_paragraph(line[2:],style='Title')
        elif line.startswith('## '):doc.add_paragraph(line[3:],style='Heading 1')
        elif line.startswith('### '):doc.add_paragraph(line[4:],style='Heading 2')
        elif line.startswith('- '):formatted(doc,line[2:],style='List Bullet')
        else:formatted(doc,line)
        if line.startswith('**Supplementary Figure'):
            key=re.match(r'\*\*Supplementary Figure (S\d+)',line).group(1);captions[key]=line
        i+=1
    primary=BASE/'Aging_Proteomics_Manuscript_Revised.docx';save_document(doc,primary,compact_images=True)
    for target in [ROOT/'Aging_Hallmarks_Proteomics_Report.docx',ROOT/'results/Aging_Hallmarks_Proteomics_Report.docx',BASE/'Aging_Proteomics_Manuscript_Final.docx']:
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(primary,target)
    supplement=styled_document();supplement.add_paragraph('Supplementary figures',style='Title')
    add_image(supplement,BASE/'figures/figS1_qc_and_null_calibration.png');formatted(supplement,captions['S1'])
    parts=sorted((BASE/'figures').glob('figS2_atlas_part_*.png'))
    for j,path in enumerate(parts,1):
        supplement.add_page_break()
        supplement.add_paragraph(f'Supplementary Figure S2 part {j} of {len(parts)}',style='Heading 1')
        add_image(supplement,path)
        if j==1:formatted(supplement,captions['S2'])
        else:formatted(supplement,'Continuation of the full alphabetical response atlas. Asterisks denote global BH q < 0.05 (*) and q < 0.01 (**); colors are clipped at +/-10.')
    supplement.add_page_break();add_image(supplement,BASE/'figures/figS3_benchmarks_and_sensitivity.png');formatted(supplement,captions['S3'])
    supplement.add_page_break();add_image(supplement,BASE/'figures/figS4_layout_context.png');formatted(supplement,captions['S4'])
    save_document(supplement,BASE/'Aging_Proteomics_Supplementary_Figures.docx')
    vignette=styled_document();vignette.add_paragraph('Supplementary Figure S5',style='Title')
    add_image(vignette,BASE/'figures/figS5_four_axis_candidates.png');formatted(vignette,captions['S5'])
    save_document(vignette,BASE/'Aging_Proteomics_Supplementary_Figure_S5.docx')
    print('Canonical manuscript, synchronized Word copies, S1-S4 supplement, and separate S5 exported.',flush=True)

if __name__=='__main__':main()
