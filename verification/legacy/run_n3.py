import time, dvg, solver
t0 = time.time()
for n in (2, 3):
    res = {}
    for arcs in dvg.all_digraphs(n):
        ss = dvg.scc_symmetric(n, arcs)
        good = [U for U in solver.all_graphs_on(n) if solver.has_potential(n, arcs, U)]
        res.setdefault((ss, bool(good)), []).append((arcs, [sorted(U) for U in good]))
    for k, v in res.items():
        print(n, 'SCC-sym' if k[0] else 'non-SCC-sym', 'potential' if k[1] else 'none', len(v))
        if k[0]:
            for arcs, Us in v: print('    ', arcs, 'U options:', Us)
print(time.time() - t0, 's')
