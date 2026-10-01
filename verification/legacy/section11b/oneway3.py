"""One portal u->w: collect outcome table over matching-number features, many graph families, starts incl. u, w."""
import random, itertools, pickle
from functools import lru_cache
import networkx as nx
from collections import defaultdict

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

def random_graph(rng):
    n = rng.randint(3, 11)
    kind = rng.choice(['gnp', 'tree', 'bip', 'pend', 'dense'])
    if kind == 'gnp':
        G = nx.gnp_random_graph(n, rng.uniform(0.15, 0.6), seed=rng.randint(0, 10**9))
    elif kind == 'dense':
        G = nx.gnp_random_graph(n, rng.uniform(0.6, 0.9), seed=rng.randint(0, 10**9))
    elif kind == 'tree':
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 10**9)) if hasattr(nx, 'random_labeled_tree') else nx.random_tree(n, seed=rng.randint(0, 10**9))
        for _ in range(rng.randint(0, 2)):
            a, b = rng.sample(range(n), 2); G.add_edge(a, b)
    elif kind == 'bip':
        a = rng.randint(1, n - 1)
        G = nx.bipartite.random_graph(a, n - a, rng.uniform(0.2, 0.7), seed=rng.randint(0, 10**9))
        G = nx.convert_node_labels_to_integers(G)
    else:
        G = nx.gnp_random_graph(n, rng.uniform(0.2, 0.5), seed=rng.randint(0, 10**9))
        m = n
        for v in list(G.nodes()):
            if rng.random() < 0.3 and m < 13:
                G.add_edge(v, m); m += 1
    return G

def features(G, Gp, V, s, u, w):
    base = nu(G, V)
    sets = [(s,), (u,), (w,), (s, u), (s, w), (u, w), (s, u, w)]
    f1 = tuple(nu(G, V - set(X)) - base for X in sets)
    f2 = tuple(nu(Gp, V - set(X)) - base for X in [(), (s,), (s, u), (s, w)])
    return f1 + f2

if __name__ == '__main__':
    rng = random.Random(2024)
    table = defaultdict(set); ex = {}
    cnt = 0
    for trial in range(12000):
        G = random_graph(rng)
        n = G.number_of_nodes()
        if n < 3 or not nx.is_connected(G): continue
        nonedges = [(a, b) for a in range(n) for b in range(n) if a != b and not G.has_edge(a, b)]
        if not nonedges: continue
        u, w = rng.choice(nonedges)
        out = {v: set(G.neighbors(v)) for v in range(n)}
        out[u].add(w)
        win = dvg_win_factory({v: frozenset(x) for v, x in out.items()})
        V = frozenset(range(n))
        Gp = G.copy(); Gp.add_edge(u, w)
        for s in range(n):
            typ = 'u' if s == u else ('w' if s == w else 'o')
            d = win(V, s)
            key = (typ,) + features(G, Gp, V, s, u, w)
            table[key].add(d); cnt += 1
            ex.setdefault((key, d), (sorted(G.edges()), u, w, s))
    conf = [k for k, v in table.items() if len(v) > 1]
    print('instances', cnt, 'classes', len(table), 'conflicts', len(conf))
    for k in conf[:10]:
        print('CONFLICT', k, ex.get((k, True)), ex.get((k, False)))
    pickle.dump((dict(table), ex), open('oneway_table.pkl', 'wb'))
    for k in sorted(table):
        print(k[0], k[1:8], '|', k[8:], '->', table[k])
