import sympy as sp
from sympy import QQ, field
K_, z = field('z', QQ)
def kappa(sizes_even, sizes_odd, m, j=0):
    sizes = [sizes_even if i % 2 == 0 else sizes_odd for i in range(m)]
    if sizes[j] == 0: return K_(1)
    b = []; k = 0
    while True:
        bk = sizes[(j + k + 1) % m] - (k + 1) // m
        b.append(max(bk, 0))
        if bk <= 0: break
        k += 1
    L = len(b); Kc = [None] * (L + 2); Kc[L] = K_(1); Kc[L + 1] = K_(0)
    for k in range(L - 1, -1, -1):
        Kc[k] = z * Kc[k + 1] - b[k] * Kc[k + 2]
    return Kc[0]
for m in (4, 6):
    print('m =', m)
    for s in range(0, 4):
        for t in range(0, 4):
            k0 = kappa(s, t, m); k1 = kappa(t, s, m)
            f = sp.factor(sp.sympify((k0 / k1).as_expr()))
            print('  s=%d t=%d  kappa0 = %s   f = kappa0(s,t)/kappa0(t,s) = %s' % (s, t, sp.factor(sp.sympify(k0.as_expr())), f))
