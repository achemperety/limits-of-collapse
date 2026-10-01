"""Proposition 10.10: the directed triangle with two frozen vertices and the smallest graph U (the six edges
between C = {0, 1, 2} and F = {3, 4}). Recomputes, over Q(z) with saturation, the minimal polynomial of each
environment weight, compares it with the sextic of Proposition 10.10, and checks the sum-of-squares form
    (z^2 - 4) a^6 + 12 a^4 + 4 z a^3 - 12 a^2 + 8 = (z a^3 + 2)^2 + 4 (1 - a^2)^3,
which has no real root a for real |z| >= 4.  Sturm counts at a few rational z are printed as a cross-check.
"""
from sage.all__sagemath_singular import PolynomialRing, QQ
import dvg

n, nf = 3, 2
arcs = [(0, 1), (1, 2), (2, 0)]
N = n + nf
F = frozenset(range(n, N))
P = dvg.reachable_positions(n, arcs)
R = dvg.make_resolvent(n, arcs)
U = frozenset([(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)])

A = PolynomialRing(dvg.Kz, ['a%d' % i for i in range(N)], order='lex')
al = A.gens()
memo = {}
eqs = [dvg.mu_weighted((W - {v}) | F, U, al, A, memo) - A(R(W, v)) * dvg.mu_weighted(W | F, U, al, A, memo)
       for (W, v) in P]
I = A.ideal(eqs)
for W in {W for (W, v) in P}:
    I = I.saturation(A.ideal([dvg.mu_weighted(W | F, U, al, A, memo)]))[0]
print('solution set dimension:', I.dimension())

Ka = PolynomialRing(dvg.Kz, 'a')
a = Ka.gen()
z = dvg.Kz(dvg.z)
S = (z**2 - 4) * a**6 + 12 * a**4 + 4 * z * a**3 - 12 * a**2 + 8
for var in range(N):
    g = I.elimination_ideal([al[j] for j in range(N) if j != var]).gens()[0]
    gu = Ka(g.univariate_polynomial().list())
    tag = 'environment' if var >= n else 'cycle'
    line = '  weight a%d (%s): minimal polynomial of degree %d' % (var, tag, gu.degree())
    if var >= n:
        line += ' | equals the sextic up to a unit: %s' % (gu.monic() == S.monic())
    print(line)
print('sum-of-squares identity holds:', S == (z * a**3 + 2)**2 + 4 * (1 - a**2)**3)

try:
    import sympy
    x = sympy.symbols('x')
    for zv in (4, 5, 10, 100, 1000):
        f = sympy.Poly((zv**2 - 4) * x**6 + 12 * x**4 + 4 * zv * x**3 - 12 * x**2 + 8, x)
        print('  z = %4d: real roots (Sturm count) = %d' % (zv, f.count_roots()))
except ImportError:
    print('sympy not available; skipping the Sturm cross-check')
