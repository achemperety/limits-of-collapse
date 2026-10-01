"""Blow-up C_m(s_0..s_{m-1}) with zero self-energies: set potential <=> for every c <= min(s),
Gamma_j(cm) = K^{(j)}_{cm} / K^{(j)}_0 is independent of the start class j, where K^{(j)} is the backward
continuant of the forced play with branching numbers b_k = s_{j+k+1} - floor((k+1)/m)."""
import sympy as sp, itertools, sys
z = sp.Symbol('z')
def branching(sizes, j):
    m = len(sizes); b = []
    k = 0
    while True:
        bk = sizes[(j + k + 1) % m] - (k + 1) // m
        b.append(bk)
        if bk == 0: break
        k += 1
    return b
def backward_continuants(b):
    L = len(b)             # levels 0..L-1; b[L-1] = 0 ends the game
    K = [None] * (L + 2); K[L] = sp.Integer(1); K[L + 1] = sp.Integer(0)
    for k in range(L - 1, -1, -1):
        K[k] = sp.expand(z * K[k + 1] - b[k] * K[k + 2])
    return K
def set_potential(sizes):
    m = len(sizes); cmax = min(sizes)
    Ks = [backward_continuants(branching(sizes, j)) for j in range(m)]
    for c in range(1, cmax + 1):
        vals = [sp.cancel(K[c * m] / K[0]) for K in Ks]
        if any(sp.simplify(v - vals[0]) != 0 for v in vals[1:]):
            return False
    return True
tests = [(1,2,1,2), (1,3,1,3), (2,3,2,3), (1,2,1,2,1,2), (2,1,1,2,1,1), (1,2,1,1,2,1), (3,1,3,1), (1,1,3,3),
         (2,5,2,5), (3,4,3,4,3,4), (1,4,1,4,1,4,1,4), (2,2,2), (2,3,2), (1,2,3,1,2,3), (4,1,4,1,4,1)]
for t in tests:
    print(t, set_potential(list(t)))
