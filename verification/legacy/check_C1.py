"""V1: in every potential found at n<=4 (over the algebraic closure of Q(z)), is every U-edge either a symmetric pair
inside one strong component or a pair of incomparable vertices?  (Theorem C proves this for real weights.)"""
import pickle, dvg
details = pickle.load(open('n4_details.pkl', 'rb'))
tot = viol = 0
kinds = {}
for arcs, ss, good in details:
    A = set(arcs); reach = dvg.reach_matrix(4, arcs)
    for U in good:
        tot += 1
        for a, b in U:
            if (a, b) in A and (b, a) in A and reach[a][b] and reach[b][a]: k = 'sym-in-comp'
            elif not reach[a][b] and not reach[b][a]: k = 'incomparable'
            else: k = 'OTHER'; viol += 1
            kinds[k] = kinds.get(k, 0) + 1
print('potentials (D,U) found at n=4:', tot, '| U-edge kinds:', kinds, '| violations:', viol)
