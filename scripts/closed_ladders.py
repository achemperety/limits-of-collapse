"""Section 14.3: reachable positions, canonical subgames and distinct exact resolvents on closed ladders."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from dvgcollapse.treewidth import canonical_subgame, closed_ladder

M = int(sys.argv[1]) if len(sys.argv) > 1 else 7
for m in range(2, M + 1):
    t0 = time.time()
    D = closed_ladder(m)
    R = D.resolvent_function()
    P = D.reachable_positions()
    subgames = {canonical_subgame(D, W, v) for (W, v) in P}
    vals = {R(W, v) for (W, v) in P}
    pred = (8 * m ** 3 - 30 * m ** 2 + 94 * m - 84) // 3
    print('m=%d: positions %6d, canonical subgames %5d (cubic fit %5d), distinct resolvents %5d  %.1fs'
          % (m, len(P), len(subgames), pred, len(vals), time.time() - t0), flush=True)
