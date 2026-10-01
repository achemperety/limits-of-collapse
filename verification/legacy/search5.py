"""n=5 search for local weighted matching potentials on non-SCC-symmetric digraphs.
Pruning (proved lemmas, real weights): U has no one-way arcs of D (sum-of-squares lemma), and U contains
at least one edge between non-adjacent vertices (otherwise all weights are normal and balance forces SCC-symmetry).
Fast filter: a vertex without phantom neighbours has weight z - sum of exit resolvents, the same at every
reachable position (checked exactly at z0). Survivors: rigorous Groebner/saturation over Q(z)."""
import pickle, itertools, sys, time
from fractions import Fraction
import dvg, solver

digs = pickle.load(open('digraphs5.pkl', 'rb'))
n = 5
part, nparts = int(sys.argv[1]), int(sys.argv[2])
z0 = Fraction(37, 11)

def resolvent_num(out, z0):
    memo = {}
    def R(W, v):
        k = (W, v)
        if k not in memo:
            s = Fraction(0)
            for w in out[v]:
                if w in W and w != v: s += R(W - {v}, w)
            memo[k] = 1 / (z0 - s)
        return memo[k]
    return R

stats = {'nonsym': 0, 'cands': 0, 'fast_pass': 0, 'groebner_pass': 0}
found = []
t0 = time.time()
for di, arcs in enumerate(digs):
    if di % nparts != part: continue
    if dvg.scc_symmetric(n, arcs): continue
    stats['nonsym'] += 1
    A = set(arcs)
    out = dvg.out_nbrs(n, arcs)
    sym = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) in A and (b, a) in A]
    nonadj = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) not in A and (b, a) not in A]
    if not nonadj: continue
    P = dvg.reachable_positions(n, arcs)
    Rn = resolvent_num(out, z0)
    byv = {v: [W for (W, u) in P if u == v] for v in range(n)}
    for ks in range(len(sym) + 1):
        for Us in itertools.combinations(sym, ks):
            for kn in range(1, len(nonadj) + 1):
                for Un in itertools.combinations(nonadj, kn):
                    U = frozenset(Us + Un)
                    stats['cands'] += 1
                    NU = {v: set() for v in range(n)}
                    for a, b in U: NU[a].add(b); NU[b].add(a)
                    ok = True
                    for v in range(n):
                        if NU[v] - out[v]: continue  # has phantom neighbours
                        exits = out[v] - NU[v]
                        vals = set()
                        for W in byv[v]:
                            vals.add(z0 - sum((Rn(W - {v}, t) for t in exits if t in W), Fraction(0)))
                            if len(vals) > 1: ok = False; break
                        if not ok: break
                    if not ok: continue
                    stats['fast_pass'] += 1
                    if solver.has_potential(n, arcs, U):
                        stats['groebner_pass'] += 1
                        found.append((arcs, sorted(U)))
                        print('FOUND', arcs, sorted(U), flush=True)
    if stats['nonsym'] % 200 == 0:
        print(part, di, stats, round(time.time() - t0), flush=True)
print('DONE', part, stats, round(time.time() - t0), flush=True)
pickle.dump(found, open('found5_%d.pkl' % part, 'wb'))
