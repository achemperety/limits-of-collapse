import itertools
from functools import lru_cache
import networkx as nx
E = [(0,8),(0,9),(1,8),(1,9),(2,6),(2,7),(2,8),(3,7),(3,8),(4,6),(4,8),(5,6),(5,7),(5,9)]
arcs = [(1,4),(2,0)]
G = nx.Graph(); G.add_nodes_from(range(10)); G.add_edges_from(E)
out = {v: frozenset(set(G.neighbors(v)) | {b for a, b in arcs if a == v}) for v in range(10)}
@lru_cache(maxsize=None)
def win(W, v): return any(not win(W - {v}, t) for t in out[v] if t in W and t != v)
nu = lambda S: len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True))
V = frozenset(range(10)); base = nu(V)
def feats(s):
    pts = [s, 1, 4, 2, 0]
    return tuple(nu(V - {pts[i] for i in X}) - base for r in range(1, 6) for X in itertools.combinations(range(5), r))
print('nu(G) =', base)
print('start 5: outcome', win(V, 5), '| start 3: outcome', win(V, 3))
print('identical deficit vectors:', feats(5) == feats(3))
print('deficits for start 5:', feats(5))
