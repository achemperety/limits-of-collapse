"""Markdown version of the referee report: resolve references and custom macros, then pandoc."""
import json, re, subprocess
refs = json.load(open('paper_refs.json'))
s = open('referee_report.tex').read()
s = re.sub(r'\\pref\{([^}]*)\}', lambda m: '%s %s' % tuple(refs[m.group(1)]), s)
s = re.sub(r'\\pnum\{([^}]*)\}', lambda m: refs[m.group(1)][1], s)
# explicit item labels
for part in ('I', 'II', 'III'):
    pat = r'\\begin\{enumerate\}\[label=\\textbf\{%s\.\\arabic\*\}[^\]]*\]' % part
    m = re.search(pat, s)
    if not m:
        continue
    start = m.end()
    end = s.index(r'\end{enumerate}', start)
    body = s[start:end]
    # only top-level \item (not inside the nested tabular)
    k = [0]
    def lab(mm):
        k[0] += 1
        return r'\item \textbf{%s.%d} ' % (part, k[0])
    body = re.sub(r'\\item ', lab, body)
    s = s[:m.start()] + r'\begin{itemize}' + body + r'\end{itemize}' + s[end + len(r'\end{enumerate}'):]
macros = {r'\Qzbar': r'\overline{\mathbb{Q}(z)}', r'\Qz': r'\mathbb{Q}(z)', r'\Rz': r'\mathbb{R}(z)',
          r'\Puis': r'\mathcal{P}_{\mathbb{R}}', r'\PSPACE': r'\mathsf{PSPACE}', r'\Pclass': r'\mathsf{P}',
          r'\poly': r'\operatorname{poly}'}
for k, v in macros.items():
    s = s.replace(k + '}', v + '}').replace(k + '$', v + '$').replace(k + ' ', v + ' ').replace(k + '(', v + '(').replace(k + ',', v + ',').replace(k + ')', v + ')').replace(k + '^', v + '^').replace(k + '=', v + '=')
s = re.sub(r'\\sev\{([^}]*)\}', r'\\textbf{[\1]}', s)
s = s.replace(r'\fix', r'\par\noindent\emph{Resolution.} ')
s = re.sub(r'\\where\{', r'\\par\\noindent\\emph{Location.} {', s)
body = s[s.index(r'\begin{document}') + len(r'\begin{document}'):s.index(r'\end{document}')]
body = body.replace(r'\maketitle', '')
body = re.sub(r'\\setlist\{[^}]*\}', '', body)
open('_report_md.tex', 'w').write('\\documentclass{article}\n\\begin{document}\n' + body + '\n\\end{document}\n')
subprocess.run(['pandoc', '_report_md.tex', '-f', 'latex', '-t', 'gfm', '--wrap=none', '-o', 'referee_report.md'], check=True)
md = open('referee_report.md').read()
md = '# Referee report (Phase I) and response map\n\n*September 30, 2026. Numbers refer to the revised paper (siamart version).*\n\n' + md
open('referee_report.md', 'w').write(md)
print(len(md.split()), 'words')
