import time, sys, dvg, solver, pickle
t0 = time.time()
n = 4
res = {}; details = []
for idx, arcs in enumerate(dvg.all_digraphs(n)):
    ss = dvg.scc_symmetric(n, arcs)
    good = []
    for U in solver.all_graphs_on(n):
        if solver.has_potential(n, arcs, U):
            good.append(sorted(U))
    res.setdefault((ss, bool(good)), 0); res[(ss, bool(good))] += 1
    details.append((arcs, ss, good))
    if idx % 20 == 0: print(idx, time.time() - t0, flush=True)
print(res, time.time() - t0, 's')
pickle.dump(details, open('n4_details.pkl', 'wb'))
