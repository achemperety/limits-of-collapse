"""Proposition 14.5(2): the comb Q_m with prime tails has 2^m distinct resolvents at the collector."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from dvgcollapse.treewidth import comb, resolvents_at_collector

M = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for m in range(1, M + 1):
    t0 = time.time()
    n = comb(m)[0].n
    k = len(resolvents_at_collector(m))
    print('m=%d n=%d distinct resolvents at q: %d (2^m = %d)  %.1fs' % (m, n, k, 2 ** m, time.time() - t0), flush=True)
