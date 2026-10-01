"""Dichotomy for two one-way arcs: the three non-star geometries admit pairs of positions with identical key-local
matching data (all nu(K - X), X a set of key points, and the adjacency among key points) but opposite outcomes."""
import itertools
from functools import lru_cache
import networkx as nx

def key_data(edges, n, arcs, s):
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(edges)
    nu = lambda S: len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True))
    V = frozenset(range(n)); base = nu(V)
    key = [s] + [x for a in arcs for x in a]
    feats = tuple(nu(V - {key[i] for i in X}) - base for r in range(1, 6) for X in itertools.combinations(range(5), r))
    pattern = tuple(tuple(i for i, x in enumerate(key) if x == y) for y in key)
    adj = tuple(G.has_edge(key[i], key[j]) for i in range(5) for j in range(i + 1, 5))
    out = {v: frozenset(set(G.neighbors(v)) | {b for a, b in arcs if a == v}) for v in range(n)}
    @lru_cache(maxsize=None)
    def win(W, v): return any(t in W and t != v and not win(W - {v}, t) for t in out[v])
    return (pattern, feats, adj), win(V, s)

pairs = {
 'chained (u1 -> w1 = u2 -> w2)': (
   ([(0, 2), (1, 2), (1, 4), (2, 5), (3, 4), (4, 5)], 6, [(0, 3), (3, 1)], 5),
   ([(0, 6), (0, 7), (0, 8), (0, 9), (1, 6), (1, 7), (1, 9), (2, 6), (2, 8), (2, 9), (3, 6), (3, 8), (3, 9), (4, 6), (4, 7),
     (4, 9), (5, 6), (5, 7), (5, 9)], 10, [(1, 3), (3, 4)], 0)),
 'common head': (
   ([(0, 5), (0, 6), (1, 5), (1, 6), (2, 6), (2, 7), (3, 5), (4, 6), (4, 7)], 8, [(1, 4), (0, 4)], 3),
   ([(0, 5), (0, 6), (0, 7), (1, 5), (1, 6), (1, 7), (2, 5), (2, 7), (3, 5), (3, 6), (4, 5), (4, 6), (4, 7)], 8, [(4, 0), (1, 0)], 2)),
 'disjoint (Proposition 10.12)': (
   ([(0, 8), (0, 9), (1, 8), (1, 9), (2, 6), (2, 7), (2, 8), (3, 7), (3, 8), (4, 6), (4, 8), (5, 6), (5, 7), (5, 9)], 10, [(1, 4), (2, 0)], 5),
   ([(0, 8), (0, 9), (1, 8), (1, 9), (2, 6), (2, 7), (2, 8), (3, 7), (3, 8), (4, 6), (4, 8), (5, 6), (5, 7), (5, 9)], 10, [(1, 4), (2, 0)], 3)),
}
for name, (A, B) in pairs.items():
    ka, wa = key_data(*A); kb, wb = key_data(*B)
    print('%-32s same key data: %s | outcomes: %s vs %s' % (name, ka == kb, 'win' if wa else 'loss', 'win' if wb else 'loss'))
