"""Seven-vertex strongly connected digraphs passing the back-arc test: which have a set potential with no
self-energies, and what are their pruned cores?"""
import sys, time
from collections import Counter
from setpot_symbolic import set_potential_exact
from describe_cores import parse, symmetric_branches, blowup
from census import prune
fname = sys.argv[1]
rows = [parse(L) for L in open(fname) if L.strip()]
t0 = time.time(); found = []; kinds = Counter(); cores_pass = Counter()
for n, arcs in rows:
    if set_potential_exact(n, arcs):
        found.append((n, arcs))
        k, sub = prune(set(range(n)), arcs)
        b = blowup(k, sub)
        key = (k, str(min(b[i:] + b[:i] for i in range(len(b)))) if b else 'not a blow-up: ' + str(sub))
        kinds[key] += 1
print('strongly connected survivors:', len(rows), '| with a set potential (eta = 0):', len(found), '| %d s' % (time.time() - t0))
for key, v in sorted(kinds.items()):
    print('   pruned core size %d, %s: %d' % (key[0], key[1], v))
pruned = [(n, a) for n, a in rows if not symmetric_branches(n, a)]
ok = [(n, a) for n, a in pruned if set_potential_exact(n, a)]
print('pruned cores on %d vertices: %d | with a set potential for eta = 0: %d' % (rows[0][0], len(pruned), len(ok)))
for n, a in ok:
    print('   ', blowup(n, a), a)
