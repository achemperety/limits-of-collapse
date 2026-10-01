"""Independent check of the Gallai-Edmonds conflicts found by ge_barrier.py.
Maximum matchings by brute force (no networkx), outcomes by exhaustive search."""
import itertools
from functools import lru_cache

def all_matchings(V, E):
    E = [e for e in E if e[0] in V and e[1] in V]
    best = [0]; out = []
    def rec(i, used, cur):
        if i == len(E):
            out.append(frozenset(cur)); return
        rec(i + 1, used, cur)
        a, b = E[i]
        if a not in used and b not in used:
            rec(i + 1, used | {a, b}, cur + [E[i]])
    rec(0, frozenset(), [])
    return out

def nu_and_ge(V, E):
    ms = all_matchings(V, E)
    nu = max(len(m) for m in ms)
    maxm = [m for m in ms if len(m) == nu]
    covered = [set(x for e in m for x in e) for m in maxm]
    D = {v for v in V if any(v not in c for c in covered)}
    adj = {v: set() for v in V}
    for a, b in E:
        if a in V and b in V:
            adj[a].add(b); adj[b].add(a)
    A = {v for v in V if v not in D and adj[v] & D}
    return nu, {v: ('D' if v in D else 'A' if v in A else 'C') for v in V}

def data(n, E, arcs, s):
    V = frozenset(range(n))
    key = [s] + [x for a in arcs for x in a]
    pattern = tuple(tuple(i for i, x in enumerate(key) if x == y) for y in key)
    Es = set(map(frozenset, E))
    adj = tuple(frozenset((key[i], key[j])) in Es for i in range(len(key)) for j in range(i + 1, len(key)))
    base, _ = nu_and_ge(V, E)
    feats = []
    for r in range(len(key) + 1):
        for X in itertools.combinations(range(len(key)), r):
            rem = {key[i] for i in X}
            if len(rem) != r:
                feats.append(None); continue
            nu, cl = nu_and_ge(V - rem, E)
            feats.append((nu - base, tuple(cl.get(key[i], '-') for i in range(len(key)))))
    return pattern, adj, tuple(feats)

def outcome(n, E, arcs, s):
    out = {v: set() for v in range(n)}
    for a, b in E:
        out[a].add(b); out[b].add(a)
    for a, b in arcs:
        out[a].add(b)
    @lru_cache(maxsize=None)
    def win(W, v):
        return any(t in W and t != v and not win(W - {v}, t) for t in out[v])
    return win(frozenset(range(n)), s)

pairs = [
 ((7, [(0,1),(0,3),(0,4),(0,5),(1,4),(1,6),(2,4),(2,6),(3,4),(3,6),(4,6),(5,6)], [(5,2),(3,2)], 0),
  (5, [(0,1),(0,2),(0,3),(0,4),(1,3),(2,3)], [(1,4),(2,4)], 3)),
 ((6, [(0,1),(1,2),(1,4),(1,5),(2,3),(2,4),(2,5)], [(5,0),(4,0)], 3),
  (6, [(0,1),(0,2),(1,2),(1,3),(1,4),(2,4),(2,5)], [(3,4),(5,4)], 0)),
 ((6, [(0,1),(0,2),(0,3),(0,4),(1,3),(2,3),(3,5)], [(2,1),(1,4)], 5),
  (6, [(0,1),(0,4),(1,2),(1,3),(1,4),(2,4),(4,5)], [(3,0),(0,2)], 5)),
 ((7, [(0,2),(0,3),(0,4),(0,6),(1,2),(1,3),(1,5),(3,4),(3,5)], [(6,2),(2,4)], 1),
  (7, [(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,4),(2,6),(3,4),(3,6),(4,5),(4,6)], [(5,2),(2,1)], 6)),
]
for P, Q in pairs:
    for (n, E, arcs, s) in (P, Q):
        Es = set(map(frozenset, E))
        assert all(frozenset(a) not in Es for a in arcs)
    dP, dQ = data(*P), data(*Q)
    print('same GE-local data:', dP == dQ, '| outcomes:', outcome(*P), outcome(*Q))
