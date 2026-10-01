"""Exact check of the two path-data identities used in the proof, on concrete digraphs.
For a one-way arc x->y and a directed path y=q0->...->qk=x, with R={q1..q_{k-1}}:
  kappa(empty) = lam(x)+lam(y)-lam(xy)          should be 0
  kappa(R)     = lam(R+x)+lam(R+y)-lam(T)-lam(R) should be 1
where gamma_X = z^|X| * Gamma(P) = 1 + lam_X z^-2 + O(z^-3) for a directed path P with V(P)=X."""
import random, itertools, sympy as sp
from functools import lru_cache
z = sp.Symbol('z')

def game(n, arcs):
    out = {v: [b for a, b in arcs if a == v] for v in range(n)}
    @lru_cache(maxsize=None)
    def R(W, v):
        s = sum((R(W - {v}, t) for t in out[v] if t in W and t != v), sp.Integer(0))
        return sp.cancel(1 / (z - s))
    def Gamma(path):
        W = frozenset(range(n)); g = sp.Integer(1)
        for t in path:
            g *= R(W, t); W = W - {t}
        return sp.cancel(g)
    return Gamma

def lam(Gamma, path):
    g = sp.cancel(z ** len(path) * Gamma(path))
    s = sp.Symbol('s')
    ser = sp.series(g.subs(z, 1 / s), s, 0, 4).removeO()
    c0, c1, c2 = [ser.coeff(s, k) for k in range(3)]
    assert c0 == 1 and c1 == 0, (c0, c1)
    return c2

def check(n, arcs, x, y, qpath):
    G = game(n, arcs)
    k = len(qpath) - 1
    assert qpath[0] == y and qpath[-1] == x
    Rset = qpath[1:-1]
    k0 = lam(G, [x]) + lam(G, [y]) - lam(G, [x, y])
    k1 = lam(G, qpath[1:]) + lam(G, qpath[:-1]) - lam(G, qpath) - (lam(G, Rset) if Rset else 0)
    return k0, k1

tests = [
    ('directed triangle', 3, [(0,1),(1,2),(2,0)], 2, 0, [0,1,2]),
    ('C4 with backward chord 2->0', 4, [(0,1),(1,2),(2,3),(3,0),(2,0)], 3, 0, [0,1,2,3]),
    ('C5 plus symmetric pair 1<->2', 5, [(0,1),(1,2),(2,1),(2,3),(3,4),(4,0)], 4, 0, [0,1,2,3,4]),
    ('two-vertex-exit C4', 6, [(0,1),(1,2),(2,3),(3,0),(0,4),(2,5),(4,5)], 3, 0, [0,1,2,3]),
]
for name, n, arcs, x, y, q in tests:
    print('%-32s kappa(empty)=%s  kappa(R)=%s' % (name, *check(n, arcs, x, y, q)))

# random strongly connected digraphs with a one-way arc on a cycle
rng = random.Random(7)
def scc_path(n, arcs, y, x):
    out = {v: [b for a, b in arcs if a == v] for v in range(n)}
    prev = {y: None}; st = [y]
    while st:
        u = st.pop(0)
        for w in out[u]:
            if w not in prev: prev[w] = u; st.append(w)
    if x not in prev: return None
    p = [x]
    while p[-1] != y: p.append(prev[p[-1]])
    return p[::-1]
cnt = 0; bad = 0
while cnt < 25:
    n = rng.randint(4, 6)
    arcs = sorted({(a, b) for a in range(n) for b in range(n) if a != b and rng.random() < 0.35})
    A = set(arcs)
    oneway = [(a, b) for a, b in arcs if (b, a) not in A]
    rng.shuffle(oneway)
    for (x, y) in oneway:
        q = scc_path(n, arcs, y, x)
        if q and len(q) >= 3:
            k0, k1 = check(n, arcs, x, y, q)
            cnt += 1; bad += (k0 != 0 or k1 != 1)
            break
print('random instances:', cnt, ' violations of (0,1):', bad)
