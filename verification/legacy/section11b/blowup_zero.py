"""Blow-ups of directed cycles with zero self-energies: set potential or not (exact test over Q(z))."""
import sys, time
from setpot_symbolic import set_potential_exact
def blowup_arcs(sizes):
    classes, v = [], 0
    for s in sizes:
        classes.append(list(range(v, v + s))); v += s
    m = len(sizes)
    arcs = [(a, b) for i in range(m) for a in classes[i] for b in classes[(i + 1) % m]]
    return v, arcs
tests = [(1,1,1), (2,2,2), (3,3,3), (1,1,2), (1,2,2), (2,2,3), (1,1,1,1), (1,2,1,2), (2,1,2,1), (2,2,2,2), (1,3,1,3),
         (2,3,2,3), (1,1,2,2), (1,2,2,1), (1,1,1,2), (1,1,1,1,1), (2,2,2,2,2), (1,1,1,1,2), (1,2,1,2,1,2),
         (1,1,1,1,1,1), (2,1,1,2,1,1), (1,2,1,1,2,1), (3,1,3,1), (1,1,3,3)]
for sizes in tests:
    n, arcs = blowup_arcs(sizes)
    t0 = time.time()
    if n > 10: print(sizes, 'skipped (n = %d)' % n); continue
    print(sizes, 'n=%d' % n, 'set potential with zero self-energies:', set_potential_exact(n, arcs), '(%.1f s)' % (time.time() - t0), flush=True)
