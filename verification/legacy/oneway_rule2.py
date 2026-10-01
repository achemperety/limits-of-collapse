"""Verify the simplified one-portal theorem on every reachable position (random graphs, pendant-leaf graphs,
trees, bipartite graphs, dense graphs). Rule:
  s = w : w essential in H
  s = u : u essential in H + uw
  s other, s essential in H : nu(H-s-u-w) < nu(H-s)  or  u inessential in H
  s other, s inessential    : nu(H-u-w) = nu(H)  and  nu(H-s-u) < nu(H)
"""
import random, sys
from functools import lru_cache
import networkx as nx
from oneway3 import random_graph

def nu_factory(G):
    @lru_cache(maxsize=None)
    def nu(S):
        return len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True))
    return nu

def rule(nu, W, s, u, w):
    b = nu(W)
    ess = lambda x: nu(W - {x}) < b
    if s == w: return ess(w)
    if s == u:  # u essential in H+uw:  nu(H-u) < nu(H+uw)
        return nu(W - {u}) < b + (1 if nu(W - {u, w}) == b else 0)
    if ess(s): return nu(W - {s, u, w}) < nu(W - {s}) or not ess(u)
    return nu(W - {u, w}) == b and nu(W - {s, u}) < b

def dvg_win_factory(out):
    @lru_cache(maxsize=None)
    def win(W, v):
        for t in out[v]:
            if t in W and t != v and not win(W - {v}, t):
                return True
        return False
    return win

seed = int(sys.argv[1]); trials = int(sys.argv[2]); nmax = int(sys.argv[3])
rng = random.Random(seed)
checked = bad = 0
for trial in range(trials):
    G = random_graph(rng)
    n = G.number_of_nodes()
    if n < 3 or n > nmax: continue
    nonedges = [(a, b) for a in range(n) for b in range(n) if a != b and not G.has_edge(a, b)]
    if not nonedges: continue
    u, w = rng.choice(nonedges)
    out = {v: frozenset(set(G.neighbors(v)) | ({w} if v == u else set())) for v in range(n)}
    win = dvg_win_factory(out); nu = nu_factory(G)
    V = frozenset(range(n))
    seen = set(); stack = [(V, s) for s in range(n)]
    while stack:
        W, v = stack.pop()
        if (W, v) in seen: continue
        seen.add((W, v))
        for t in out[v]:
            if t in W and t != v: stack.append((W - {v}, t))
    for (W, v) in seen:
        if u in W and w in W:
            checked += 1
            if win(W, v) != rule(nu, W, v, u, w):
                bad += 1
                if bad <= 5: print('MISMATCH', sorted(G.edges()), (u, w), sorted(W), v)
print('seed', seed, 'positions checked', checked, 'mismatches', bad)
