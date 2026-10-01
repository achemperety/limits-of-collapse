import dvg
import pickle
from collections import Counter
rows = pickle.load(open('setpot_rows.pkl', 'rb'))
cnt = Counter()
for n, arcs, comps, info in rows:
    if n != 5: continue
    d = info[0]
    others = 5 - d['size']
    if d['cycle']:
        seq = [str(e) for e in d['eta_seq']]
        kind = 'C%d' % d['size']
        eta = seq[0] if len(set(seq)) == 1 else '(' + ', '.join(seq) + ')'
    else:
        kind = {((0, 1), (0, 2), (1, 0), (2, 3), (3, 0)): 'triangle + pendant',
                ((0, 1), (0, 2), (1, 3), (2, 3), (3, 0)): '2-page book',
                ((0, 1), (0, 2), (1, 0), (2, 3), (2, 4), (3, 0), (4, 0)): '2-page book + pendant at spine tail'}[d['canon']]
        eta = str({k: str(v) for k, v in d['eta'].items()})
    cnt[(kind, eta)] += 1
for k, v in sorted(cnt.items()): print(v, k)
print(sum(cnt.values()))
