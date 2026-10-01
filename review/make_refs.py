"""Resolve references to the revised paper from its .aux file (siamart build) into plain text macros."""
import re, sys
aux = open(sys.argv[1] if len(sys.argv) > 1 else 'paper.aux').read()
names = {'theorem': 'Theorem', 'lemma': 'Lemma', 'corollary': 'Corollary', 'proposition': 'Proposition',
         'definition': 'Definition', 'conjecture': 'Conjecture', 'problem': 'Open Problem', 'remark': 'Remark',
         'example': 'Example', 'claim': 'Claim', 'section': 'Section', 'subsection': 'Section',
         'appendix': 'Appendix', 'equation': 'Equation', 'figure': 'Figure'}
refs = {}
for m in re.finditer(r'\\newlabel\{([^}@]+)@cref\}\{\{\[([^\]]*)\]\[[^\]]*\]\[[^\]]*\]([^}]*)\}', aux):
    lab, typ, num = m.group(1), m.group(2), m.group(3)
    refs[lab] = (names.get(typ, typ.capitalize()), num)
out = []
for lab, (typ, num) in sorted(refs.items()):
    out.append('\\expandafter\\def\\csname pref@%s\\endcsname{%s~%s}' % (lab, typ, num))
    out.append('\\expandafter\\def\\csname pnum@%s\\endcsname{%s}' % (lab, num))
open('paper_refs.tex', 'w').write('\n'.join(out) + '\n')
import json
json.dump(refs, open('paper_refs.json', 'w'), indent=0)
print(len(refs), 'labels')
