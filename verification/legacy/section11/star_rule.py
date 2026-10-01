"""Conjecture: one-way arcs u -> w_1, ..., u -> w_k with a common tail (a 'star') obey the rule of Theorem 10.11
with the single head w replaced by the set of live heads. Checked on every reachable position with u unvisited."""
import random, sys
from functools import lru_cache
import networkx as nx
from oneway3 import random_graph

def nu_factory(G):
    @lru_cache(maxsize=None)
    def nu(S): return len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True))
    return nu

def rule(nu, nu_plus, W, s, u, heads):
    """heads: live heads (unvisited, != s). nu_plus(S, hs) = matching number of H[S] + u*hs."""
    b = nu(W)
    if s == u:
        return nu(W - {u}) < nu_plus(W, heads)
    if not heads:                       # plain undirected geography
        return nu(W - {s}) < b
    ess_s = nu(W - {s}) < b
    if ess_s:
        split = nu_plus(W - {s}, heads) == nu(W - {s})
        return split or not (nu(W - {u}) < b)
    joined = nu_plus(W, heads) > b
    return joined and nu(W - {s, u}) < b

seed, trials, nmax, kmax = map(int, sys.argv[1:5])
rng = random.Random(seed)
checked = bad = 0
for trial in range(trials):
    G = random_graph(rng); n = G.number_of_nodes()
    if n < 3 or n > nmax: continue
    u = rng.randrange(n)
    cand = [w for w in range(n) if w != u and not G.has_edge(u, w)]
    if not cand: continue
    k = rng.randint(1, min(kmax, len(cand)))
    heads = frozenset(rng.sample(cand, k))
    out = {v: frozenset(set(G.neighbors(v)) | (heads if v == u else set())) for v in range(n)}
    @lru_cache(maxsize=None)
    def win(W, v): return any(t in W and t != v and not win(W - {v}, t) for t in out[v])
    nu = nu_factory(G)
    @lru_cache(maxsize=None)
    def nu_plus(S, hs):
        H = G.subgraph(S).copy()
        if u in S:
            for w in hs:
                if w in S: H.add_edge(u, w)
        return len(nx.max_weight_matching(H, maxcardinality=True))
    V = frozenset(range(n))
    seen = set(); stack = [(V, s) for s in range(n)]
    while stack:
        W, v = stack.pop()
        if (W, v) in seen: continue
        seen.add((W, v))
        for t in out[v]:
            if t in W and t != v: stack.append((W - {v}, t))
    for (W, v) in seen:
        if u not in W: continue
        live = frozenset(w for w in heads if w in W and w != v)
        checked += 1
        if win(W, v) != rule(nu, nu_plus, W, v, u, live):
            bad += 1
            if bad <= 5: print('MISMATCH', sorted(G.edges()), u, sorted(heads), sorted(W), v, win(W, v))
print('seed', seed, 'positions checked', checked, 'mismatches', bad)
