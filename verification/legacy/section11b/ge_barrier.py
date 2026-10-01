"""Does Gallai-Edmonds data at the key points decide two-tail geography?
GE-local data of a position: the key-local data of Section 11.1 together with, for every set X of key points and
every key point y outside X, the Gallai-Edmonds class of y in K - X (D: missed by some maximum matching,
A: neighbour of D outside D, C: the rest)."""
import itertools, random, sys
from functools import lru_cache
from collections import defaultdict, Counter
import networkx as nx

from oneway3 import random_graph


def ge_classes(G, S):
    H = G.subgraph(S)
    nu = len(nx.max_weight_matching(H, maxcardinality=True))
    D = {v for v in S if len(nx.max_weight_matching(H.subgraph(set(S) - {v}), maxcardinality=True)) == nu}
    Aset = {v for v in S if v not in D and any(u in D for u in H.neighbors(v))}
    return {v: ('D' if v in D else 'A' if v in Aset else 'C') for v in S}, nu


def data(G, arcs, s, with_ge=True):
    n = G.number_of_nodes(); V = frozenset(range(n))
    key = [s] + [x for a in arcs for x in a]
    pattern = tuple(tuple(i for i, x in enumerate(key) if x == y) for y in key)
    adj = tuple(G.has_edge(key[i], key[j]) for i in range(len(key)) for j in range(i + 1, len(key)))
    base = len(nx.max_weight_matching(G, maxcardinality=True))
    feats, ges = [], []
    for r in range(0, len(key) + 1):
        for X in itertools.combinations(range(len(key)), r):
            S = V - {key[i] for i in X}
            if r and len(S) != len(V) - len({key[i] for i in X}):
                feats.append(None); ges.append(None); continue
            cl, nu = ge_classes(G, S) if with_ge else ({}, len(nx.max_weight_matching(G.subgraph(S), maxcardinality=True)))
            feats.append(nu - base)
            ges.append(tuple(cl.get(key[i], '-') for i in range(len(key))) if with_ge else None)
    return (pattern, tuple(feats), adj, tuple(ges))


def outcome(G, arcs, s):
    n = G.number_of_nodes()
    out = {v: frozenset(set(G.neighbors(v)) | {b for a, b in arcs if a == v}) for v in range(n)}

    @lru_cache(maxsize=None)
    def win(W, v):
        return any(t in W and t != v and not win(W - {v}, t) for t in out[v])
    return win(frozenset(range(n)), s)


def geometry(arcs):
    (a1, b1), (a2, b2) = arcs
    if a1 == a2: return 'common tail'
    if b1 == b2: return 'common head'
    if b1 == a2 or b2 == a1: return 'chained'
    return 'disjoint'


if __name__ == '__main__':
    pairs = {
        'chained': (([(0, 2), (1, 2), (1, 4), (2, 5), (3, 4), (4, 5)], 6, [(0, 3), (3, 1)], 5),
                    ([(0, 6), (0, 7), (0, 8), (0, 9), (1, 6), (1, 7), (1, 9), (2, 6), (2, 8), (2, 9), (3, 6), (3, 8), (3, 9),
                      (4, 6), (4, 7), (4, 9), (5, 6), (5, 7), (5, 9)], 10, [(1, 3), (3, 4)], 0)),
        'common head': (([(0, 5), (0, 6), (1, 5), (1, 6), (2, 6), (2, 7), (3, 5), (4, 6), (4, 7)], 8, [(1, 4), (0, 4)], 3),
                        ([(0, 5), (0, 6), (0, 7), (1, 5), (1, 6), (1, 7), (2, 5), (2, 7), (3, 5), (3, 6), (4, 5), (4, 6), (4, 7)], 8,
                         [(4, 0), (1, 0)], 2)),
        'disjoint': (([(0, 8), (0, 9), (1, 8), (1, 9), (2, 6), (2, 7), (2, 8), (3, 7), (3, 8), (4, 6), (4, 8), (5, 6), (5, 7), (5, 9)],
                      10, [(1, 4), (2, 0)], 5),
                     ([(0, 8), (0, 9), (1, 8), (1, 9), (2, 6), (2, 7), (2, 8), (3, 7), (3, 8), (4, 6), (4, 8), (5, 6), (5, 7), (5, 9)],
                      10, [(1, 4), (2, 0)], 3)),
    }
    for name, (Aa, Bb) in pairs.items():
        Ga = nx.Graph(); Ga.add_nodes_from(range(Aa[1])); Ga.add_edges_from(Aa[0])
        Gb = nx.Graph(); Gb.add_nodes_from(range(Bb[1])); Gb.add_edges_from(Bb[0])
        da, db = data(Ga, Aa[2], Aa[3]), data(Gb, Bb[2], Bb[3])
        print('%-12s key-local equal: %s | GE classes equal: %s' % (name, da[:3] == db[:3], da == db))
    # sample search for conflicts under GE-local data
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    rng = random.Random(SEED)
    table = defaultdict(set); ex = {}; geo = Counter(); cnt = 0
    for trial in range(trials):
        G = random_graph(rng); n = G.number_of_nodes()
        if n < 4 or n > NMAX or not nx.is_connected(G): continue
        nonedges = [(a, b) for a in range(n) for b in range(n) if a != b and not G.has_edge(a, b)]
        if len(nonedges) < 2: continue
        arcs = rng.sample(nonedges, 2)
        if len({frozenset(a) for a in arcs}) < 2 or geometry(arcs) == 'common tail': continue
        for s in range(n):
            d = data(G, arcs, s); o = outcome(G, arcs, s); cnt += 1
            table[d].add(o); ex.setdefault((d, o), (geometry(arcs), n, sorted(G.edges()), arcs, s)); geo[geometry(arcs)] += 1
    conf = [k for k, v in table.items() if len(v) > 1]
    print('two-tail starts:', cnt, dict(geo), '| GE-local classes:', len(table), '| conflicts:', len(conf))
    for k in conf[:6]:
        print('  WIN ', ex[(k, True)]); print('  LOSS', ex[(k, False)])
