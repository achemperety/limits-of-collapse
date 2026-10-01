"""Complete n=5 check over the algebraic closure of Q(z), pruned by Theorem 10.6 (no set potential => no potential
over any field). Survivors: every graph U on 5 vertices except U inside Sym(D) (excluded by Theorem 10.1 over any
field); each (D, U) decided by a Groebner basis with saturation over Q(z)."""
import pickle, itertools, sys, time
import dvg, solver
from setpot_symbolic import set_potential_exact

digs = pickle.load(open('digraphs5.pkl', 'rb'))
n = 5
allpairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
t0 = time.time()
nonsym = 0
survivors = []
for arcs in digs:
    if dvg.scc_symmetric(n, arcs):
        continue
    nonsym += 1
    if set_potential_exact(n, arcs):
        survivors.append(arcs)
print('non-SCC-symmetric:', nonsym, '| with a set potential (exact):', len(survivors),
      round(time.time() - t0), 's', flush=True)
cands = found = 0
for idx, arcs in enumerate(survivors):
    A = set(arcs)
    sym = {(a, b) for (a, b) in allpairs if (a, b) in A and (b, a) in A}
    for mask in range(1 << len(allpairs)):
        U = frozenset(allpairs[i] for i in range(len(allpairs)) if mask >> i & 1)
        if U <= sym:
            continue
        cands += 1
        if solver.has_potential(n, arcs, U):
            found += 1
            print('FOUND', arcs, sorted(U), flush=True)
    print('digraph', idx + 1, 'of', len(survivors), 'done;', cands, 'candidates so far', round(time.time() - t0), 's',
          flush=True)
print('DONE: candidates', cands, 'potentials found', found, round(time.time() - t0), 's')
