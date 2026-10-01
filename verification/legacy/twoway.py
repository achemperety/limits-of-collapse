"""k one-way arcs (inside a connected symmetric graph): is the fresh-start outcome a function of the deficits
nu(G - X) - nu(G) over subsets X of {s} u {arc endpoints}?"""
import random, itertools, sys
from functools import lru_cache
import networkx as nx
from collections import defaultdict
from oneway3 import random_graph

def nu(G, S):
    return len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True))

def dvg_win_factory(out):
    @lru_cache(maxsize=None)
    def win(W, v):
        for t in out[v]:
            if t in W and t != v and not win(W - {v}, t):
                return True
        return False
    return win

k = int(sys.argv[1]); trials = int(sys.argv[2])
rng = random.Random(77 + k)
table = defaultdict(set); cnt = 0; ex = {}
for trial in range(trials):
    G = random_graph(rng)
    n = G.number_of_nodes()
    if n < 4 or not nx.is_connected(G): continue
    nonedges = [(a, b) for a in range(n) for b in range(n) if a != b and not G.has_edge(a, b)]
    if len(nonedges) < k: continue
    arcs = rng.sample(nonedges, k)
    # avoid both directions of the same pair
    if len({frozenset(a) for a in arcs}) < k: continue
    out = {v: set(G.neighbors(v)) for v in range(n)}
    for a, b in arcs: out[a].add(b)
    win = dvg_win_factory({v: frozenset(x) for v, x in out.items()})
    V = frozenset(range(n))
    base = nu(G, V)
    ends = [x for a in arcs for x in a]
    for s in range(n):
        key_pts = [s] + ends
        feats = []
        for r in range(1, len(key_pts) + 1):
            for X in itertools.combinations(range(len(key_pts)), r):
                S = V - {key_pts[i] for i in X}
                feats.append(nu(G, S) - base if len(S) == len(V) - len({key_pts[i] for i in X}) else None)
        # role pattern: which key points coincide
        pattern = tuple(tuple(i for i, x in enumerate(key_pts) if x == y) for y in key_pts)
        key = (pattern, tuple(feats))
        d = win(V, s)
        table[key].add(d); cnt += 1
        ex.setdefault((key, d), (sorted(G.edges()), arcs, s))
conf = [kk for kk, v in table.items() if len(v) > 1]
print('k', k, 'instances', cnt, 'classes', len(table), 'conflicts', len(conf))
for kk in conf[:5]:
    print('CONFLICT', ex.get((kk, True)), '||', ex.get((kk, False)))
