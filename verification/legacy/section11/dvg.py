"""Core utilities: digraphs, DVG resolvents over Q(z), reachable positions, weighted matching polynomials."""
import itertools
from functools import lru_cache
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField

Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Kz = FractionField(Rz)

def out_nbrs(n, arcs):
    out = {v: set() for v in range(n)}
    for a, b in arcs: out[a].add(b)
    return out

def reach_matrix(n, arcs):
    out = out_nbrs(n, arcs)
    reach = [[False]*n for _ in range(n)]
    for v in range(n):
        reach[v][v] = True; st = [v]
        while st:
            x = st.pop()
            for y in out[x]:
                if not reach[v][y]: reach[v][y] = True; st.append(y)
    return reach

def scc_symmetric(n, arcs):
    A = set(arcs); reach = reach_matrix(n, arcs)
    return all(not (reach[b][a] and (b, a) not in A) for a, b in A)

def reachable_positions(n, arcs, starts=None):
    out = out_nbrs(n, arcs); seen = set()
    def dfs(W, v):
        if (W, v) in seen: return
        seen.add((W, v))
        for w in out[v]:
            if w in W and w != v: dfs(W - {v}, w)
    for s in (range(n) if starts is None else starts):
        dfs(frozenset(range(n)), s)
    return seen

def make_resolvent(n, arcs, field=Kz, zval=None):
    out = out_nbrs(n, arcs)
    zz = field(z) if zval is None else field(zval)
    @lru_cache(maxsize=None)
    def R(W, v):
        s = field(0)
        for w in out[v]:
            if w in W and w != v: s += R(W - {v}, w)
        return 1 / (zz - s)
    return R

def mu_weighted(S, E, a, ring, memo):
    """weighted matching polynomial of U[S]; E set of (i<j) pairs; a dict/list of weights in ring"""
    S = frozenset(S)
    if S in memo: return memo[S]
    if not S: return ring(1)
    v = min(S); rest = S - {v}
    val = a[v] * mu_weighted(rest, E, a, ring, memo)
    for w in rest:
        if (min(v, w), max(v, w)) in E:
            val -= mu_weighted(rest - {w}, E, a, ring, memo)
    memo[S] = val
    return val

def canon(n, arcs):
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(sorted((p[a], p[b]) for a, b in arcs))
        if best is None or key < best: best = key
    return best

def all_digraphs(n):
    pairs = [(a, b) for a in range(n) for b in range(n) if a != b]
    seen = set()
    for mask in range(1 << len(pairs)):
        arcs = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
        c = canon(n, arcs)
        if c in seen: continue
        seen.add(c); yield list(c)
