"""Typeset the complete mathematical manuscript using Pandoc and XeLaTeX.

Requires pandoc, xelatex and the standard LaTeX packages used below.
The final PDF is written under output/pdf; build intermediates under tmp/pdfs.
"""
from pathlib import Path
import subprocess
import shutil
import re

ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / 'tmp/pdfs'
OUT = ROOT / 'output/pdf'
TMP.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
raw = (ROOT / 'manuscript.md').read_text()
title, body = raw.split('\n', 1)
title = title.lstrip('# ')
body = body.lstrip('\n')
if body.startswith('**Anonymous**'):
    body = body.split('\n', 1)[1]
def grouped_figure(match):
    path, caption = match.group(2), match.group(3)
    width = {'geometry': '0.80', 'critical_scaling': '0.72'}.get(Path(path).stem, '0.82')
    return ('\\begin{figure}[tbp]\n\\centering\n'
            '\\includegraphics[width=' + width + '\\linewidth]{' + path + '}\n'
            '\\caption{' + caption + '}\n\\end{figure}')
body = re.sub(r'!\[([^\]]*)\]\((figures/[^)]+)\)\n\n\*Figure \d+\. (.*?)\*',
              grouped_figure, body)
body = body.replace('\n## Appendix A.', '\n\\FloatBarrier\n\n## Appendix A.')
(TMP / 'body.md').write_text(body)
(TMP / 'header.tex').write_text(r'''
\usepackage{microtype}
\usepackage{graphicx}
\usepackage{float}
\providecommand{\FloatBarrier}{\clearpage}
\usepackage{fancyhdr}
\usepackage{xurl}
\usepackage{needspace}
\allowdisplaybreaks[1]
\setlength{\emergencystretch}{3em}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt plus 1pt}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Alternating projections: exact worst-case convergence}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\setlength{\headheight}{14pt}
''')
cmd = ['pandoc', str(TMP / 'body.md'), '--standalone',
       '--from=markdown+tex_math_dollars+autolink_bare_uris-implicit_figures',
       '--to=latex', '--shift-heading-level-by=-1',
       '--resource-path=' + str(ROOT),
       '--include-in-header=' + str(TMP / 'header.tex'),
       '--metadata=title:' + title,
       '--metadata=author:Anonymous',
       '--metadata=date:3 November 2026',
       '-V', 'documentclass:article', '-V', 'fontsize:11pt',
       '-V', 'papersize:a4', '-V', 'geometry:margin=25mm',
       '-V', 'mainfont=lmroman10-regular.otf',
       '-V', 'mainfontoptions=BoldFont=lmroman10-bold.otf,ItalicFont=lmroman10-italic.otf,BoldItalicFont=lmroman10-bolditalic.otf',
       '-V', 'mathfont=latinmodern-math.otf',
       '-V', 'colorlinks:true', '-V', 'urlcolor:teal',
       '-o', str(TMP / 'paper.tex')]
subprocess.run(cmd, cwd=ROOT, check=True)
tex = (TMP / 'paper.tex').read_text()
# Short tables: plain centred tabulars avoid longtable/float interaction.
tex = re.sub(r'\\begin\{longtable\}\[\]\{(@\{\}r+@\{\})\}', r'\\begin{center}\\begin{tabular}{\1}', tex)
tex = tex.replace('\\end{longtable}', '\\end{tabular}\\end{center}')
tex = tex.replace('\\endhead\n', '').replace('\\endlastfoot\n', '')
tex = tex.replace('\\bottomrule\\noalign{}\n', '')
tex = re.sub(r'(\\\\\n)(\\end\{tabular\})', r'\1\\bottomrule\n\2', tex)
(TMP / 'paper.tex').write_text(tex)
logs = []
for run in (1, 2):
    result = subprocess.run(['xelatex', '-interaction=nonstopmode',
        '-halt-on-error', '-output-directory=' + str(TMP),
        str(TMP / 'paper.tex')], cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    logs.append(f'=== XeLaTeX run {run} ===\n' + result.stdout)
    if result.returncode:
        (ROOT / 'results/pdf_build_log.txt').write_text('\n'.join(logs))
        print(result.stdout[-10000:])
        raise SystemExit(result.returncode)
(ROOT / 'results/pdf_build_log.txt').write_text('\n'.join(logs))
dest = OUT / 'exact_worst_case_alternating_projections.pdf'
shutil.copyfile(TMP / 'paper.pdf', dest)
print('Created', dest)
for line in logs[-1].splitlines():
    if any(term in line for term in ['Overfull', 'Missing character', 'Warning:', 'Output written']):
        print(line)
