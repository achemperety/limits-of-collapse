"""Reproduce the k = 2 experiment of twoway.py (seed 79, 6000 trials) and list every conflicting class with the
geometry of the two arcs (disjoint, chained, common tail, common head, antiparallel)."""
import random, itertools, sys
from functools import lru_cache
import networkx as nx
from collections import defaultdict
from oneway3 import random_graph

def nu(G, S): return len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True))
def dvg_win_factory(out):
    @lru_cache(maxsize=None)
    def win(W, v):
        return any(t in W and t != v and not win(W - {v}, t) for t in out[v])
    return win
def geometry(arcs):
    (a1, b1), (a2, b2) = arcs
    if a1 == a2: return 'common tail'
    if b1 == b2: return 'common head'
    if b1 == a2 or b2 == a1:
        return 'antiparallel' if (b1 == a2 and b2 == a1) else 'chained'
    return 'disjoint'
k, trials = 2, 6000
rng = random.Random(77 + k)
table = defaultdict(list); geo_count = defaultdict(int)
for trial in range(trials):
    G = random_graph(rng); n = G.number_of_nodes()
    if n < 4 or not nx.is_connected(G): continue
    nonedges = [(a, b) for a in range(n) for b in range(n) if a != b and not G.has_edge(a, b)]
    if len(nonedges) < k: continue
    arcs = rng.sample(nonedges, k)
    if len({frozenset(a) for a in arcs}) < k: continue
    out = {v: set(G.neighbors(v)) for v in range(n)}
    for a, b in arcs: out[a].add(b)
    win = dvg_win_factory({v: frozenset(x) for v, x in out.items()})
    V = frozenset(range(n)); base = nu(G, V); ends = [x for a in arcs for x in a]
    for s in range(n):
        key_pts = [s] + ends
        feats = []
        for r in range(1, len(key_pts) + 1):
            for X in itertools.combinations(range(len(key_pts)), r):
                S = V - {key_pts[i] for i in X}
                feats.append(nu(G, S) - base if len(S) == len(V) - len({key_pts[i] for i in X}) else None)
        pattern = tuple(tuple(i for i, x in enumerate(key_pts) if x == y) for y in key_pts)
        key = (pattern, tuple(feats))
        table[key].append((win(V, s), sorted(G.edges()), arcs, s, n))
        geo_count[geometry(arcs)] += 1
conf = {kk: v for kk, v in table.items() if len({x[0] for x in v}) > 1}
print('instances', sum(len(v) for v in table.values()), 'classes', len(table), 'conflicts', len(conf))
print('instances by arc geometry:', dict(geo_count))
for kk, v in conf.items():
    g = {geometry(x[2]) for x in v}
    starts_in_key = kk[0][0]
    print('CLASS geometry', g, '| start coincides with key points', starts_in_key, '| sizes', sorted({x[4] for x in v}),
          '| instances', len(v), 'wins', sum(x[0] for x in v))

print()
for kk, v in conf.items():
    v_sorted = sorted(v, key=lambda x: x[4])
    wins = [x for x in v_sorted if x[0]][:1]; loss = [x for x in v_sorted if not x[0]][:1]
    print('--- class', {geometry(x[2]) for x in v}, 'deficit vector', kk[1][:5], '...')
    for x in wins + loss:
        print('    %s  n=%d edges=%s arcs=%s start=%d' % ('WIN ' if x[0] else 'LOSS', x[4], x[1], x[2], x[3]))
