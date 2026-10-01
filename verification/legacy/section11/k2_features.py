"""Does adding the adjacency pattern among the key points (s and the arc endpoints) remove the k = 2 conflicts?"""
import random, itertools, sys
from functools import lru_cache
import networkx as nx
from collections import defaultdict
from oneway3 import random_graph
def geometry(arcs):
    (a1, b1), (a2, b2) = arcs
    if a1 == a2: return 'common tail'
    if b1 == b2: return 'common head'
    if b1 == a2 or b2 == a1: return 'antiparallel' if (b1 == a2 and b2 == a1) else 'chained'
    return 'disjoint'
seed, trials = int(sys.argv[1]), int(sys.argv[2])
k = 2
rng = random.Random(seed)
T0 = defaultdict(set); T1 = defaultdict(set); ex1 = {}
cnt = 0
for trial in range(trials):
    G = random_graph(rng); n = G.number_of_nodes()
    if n < 4 or not nx.is_connected(G): continue
    nonedges = [(a, b) for a in range(n) for b in range(n) if a != b and not G.has_edge(a, b)]
    if len(nonedges) < k: continue
    arcs = rng.sample(nonedges, k)
    if len({frozenset(a) for a in arcs}) < k: continue
    out = {v: set(G.neighbors(v)) for v in range(n)}
    for a, b in arcs: out[a].add(b)
    out = {v: frozenset(x) for v, x in out.items()}
    @lru_cache(maxsize=None)
    def win(W, v): return any(t in W and t != v and not win(W - {v}, t) for t in out[v])
    @lru_cache(maxsize=None)
    def nu(S): return len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True))
    V = frozenset(range(n)); base = nu(V); ends = [x for a in arcs for x in a]
    for s in range(n):
        key_pts = [s] + ends
        feats = []
        for r in range(1, len(key_pts) + 1):
            for X in itertools.combinations(range(len(key_pts)), r):
                S = V - {key_pts[i] for i in X}
                feats.append(nu(S) - base if len(S) == len(V) - len({key_pts[i] for i in X}) else None)
        pattern = tuple(tuple(i for i, x in enumerate(key_pts) if x == y) for y in key_pts)
        adj = tuple(G.has_edge(key_pts[i], key_pts[j]) for i in range(5) for j in range(i + 1, 5))
        d = win(V, s); cnt += 1
        T0[(pattern, tuple(feats))].add(d)
        key1 = (pattern, tuple(feats), adj)
        T1[key1].add(d)
        ex1.setdefault((key1, d), (geometry(arcs), n, sorted(G.edges()), arcs, s))
c0 = [x for x, v in T0.items() if len(v) > 1]; c1 = [x for x, v in T1.items() if len(v) > 1]
print('instances', cnt, '| deletion numbers only: classes', len(T0), 'conflicts', len(c0),
      '| + key adjacency: classes', len(T1), 'conflicts', len(c1))
from collections import Counter
print('conflicts after adding adjacency, by geometry:', Counter(ex1[(x, True)][0] for x in c1))
for x in c1[:4]:
    print('  WIN ', ex1[(x, True)]); print('  LOSS', ex1[(x, False)])
