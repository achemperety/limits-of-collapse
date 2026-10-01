"""Proposition 11.10 (triangle plus three isolated vertices): stage 1 filter at two rational z, stage 2 exact.
Requires passagemath-singular (Singular)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
"""Exhaustive search: D = directed triangle {0,1,2} + isolated vertices {3,4,5}; every graph U on 6 vertices up
to the automorphisms of D (rotations of the triangle x permutations of the isolated vertices).
Stage 1 (filter, rigorous for 'no potential' only generically): saturated Groebner basis with z specialised to
two rational values.  Stage 2 (rigorous): saturated Groebner basis over Q(z) for every U that survives stage 1."""
import itertools, sys, time, pickle
from dvgcollapse.dvg import Digraph
from dvgcollapse import groebner_saturation as gs
from sage.all__sagemath_singular import PolynomialRing, QQ
from fractions import Fraction

n = 6
D = Digraph(n, [(0, 1), (1, 2), (2, 0)])
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
autos = []
for r in range(3):
    for p in itertools.permutations([3, 4, 5]):
        perm = {i: (i + r) % 3 for i in range(3)}
        perm.update({3 + i: p[i] for i in range(3)})
        autos.append(perm)
def canon(U):
    best = None
    for s in autos:
        key = tuple(sorted(tuple(sorted((s[a], s[b]))) for a, b in U))
        if best is None or key < best: best = key
    return best
orbits = {}
for mask in range(1 << len(pairs)):
    U = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
    orbits.setdefault(canon(U), U)
reps = list(orbits.values())
print('graphs U up to symmetry:', len(reps), flush=True)

# quick necessary condition (any field, Theorem 10.7): the triangle lies in one component of U
def triangle_connected(U):
    comp = {v: v for v in range(n)}
    def find(x):
        while comp[x] != x: x = comp[x]
        return x
    for a, b in U: comp[find(a)] = find(b)
    return find(0) == find(1) == find(2)

# stage 1: specialised z
def decide_at(U, z0):
    pos, vals, adj, F = gs.potential_system(D, U, 0, None)
    A = PolynomialRing(QQ, ['a%d' % i for i in range(n)], order='degrevlex')
    al = A.gens(); a = {i: al[i] for i in range(n)}; one = A(1); memo = {}
    eqs = []; Ws = set()
    for (W, v), r in zip(pos, vals):
        # evaluate the rational function r at z0 exactly
        num = r.numer; den = r.denom
        def ev(p):
            s = Fraction(0)
            for (e,), c in p.terms(): s += Fraction(int(c.numerator), int(c.denominator)) * z0 ** e
            return s
        rv = ev(num) / ev(den)
        e = gs._mu((W - {v}) | F, adj, a, one, memo) - A(QQ(rv.numerator) / QQ(rv.denominator)) * gs._mu(W | F, adj, a, one, memo)
        if e != 0: eqs.append(e)
        Ws.add(W)
    nonzero = [gs._mu(W | F, adj, a, one, memo) for W in Ws]
    if any(f == 0 for f in nonzero): return False
    if not eqs: return True
    J = A.ideal(eqs)
    for f in nonzero:
        if f.is_constant(): continue
        J = J.saturation(A.ideal([f]))[0]
        if J.is_one(): return False
    return not J.is_one()

t0 = time.time(); surv = []
cand = [U for U in reps if triangle_connected(U)]
print('after Theorem 10.7 filter:', len(cand), flush=True)
for idx, U in enumerate(cand):
    if decide_at(U, Fraction(37, 11)) and decide_at(U, Fraction(-53, 7)):
        surv.append(U)
    if idx % 200 == 0: print(idx, len(surv), round(time.time() - t0), 's', flush=True)
print('stage 1 survivors:', len(surv), round(time.time() - t0), 's', flush=True)
pickle.dump(surv, open('c3k3_survivors.pkl', 'wb'))
found = []
for U in surv:
    t1 = time.time(); ok = gs.decide(D, U, backend='singular')
    print('stage 2', U, ok, round(time.time() - t1, 1), 's', flush=True)
    if ok: found.append(U)
print('RESULT: potentials on C3 + 3 isolated vertices over the algebraic closure:', len(found), found, flush=True)
