"""Sufficiency check of the single-class blow-up family for m = 3..12 with r, w, p symbolic."""
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField
def cont(seq, coup):
    a, b = 1, seq[0]
    for i in range(1, len(seq)):
        a, b = b, seq[i] * b - coup[i] * a
    return b
P = PolynomialRing(QQ, ['w', 'p', 'r']); F = FractionField(P)
w, p, r = [F(g) for g in P.gens()]
for m in range(3, 13):
    if m % 2:
        ws = [None] + [w - (r - 1) / w] * (m - 2) + [w]
    else:
        ws = [None] + [(p - (r - 1) / w) if i % 2 else (w - (r - 1) / p) for i in range(1, m - 1)] + [p]
    Ks = []
    for a in range(1, m):
        seq = [ws[i] for i in range(a, m)] + [w] + [ws[i] for i in range(1, a)]
        coup = [0] + [r if (i == m - a) else 1 for i in range(1, m)]
        Ks.append(cont(seq, coup))
    Kt = cont([w] + [ws[i] for i in range(1, m - 1)] + [ws[m - 1] - (r - 1) / w], [0] * m)
    Kt = cont([w] + [ws[i] for i in range(1, m - 1)] + [ws[m - 1] - (r - 1) / w], [0] + [1] * (m - 1))
    print('m=%2d: all rotation continuants equal: %s' % (m, all(K == Kt for K in Ks)), '| K =', Kt.numerator().factor() if m <= 5 else '')
