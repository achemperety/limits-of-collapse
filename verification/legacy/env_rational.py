"""V6: frozen-environment potentials for the directed triangle with 2 frozen vertices: are any of them rational?"""
import dvg
from sage.all__sagemath_singular import PolynomialRing, QQ
n, nf = 3, 2
arcs = [(0,1),(1,2),(2,0)]
N = n + nf
F = frozenset(range(n, N))
P = dvg.reachable_positions(n, arcs)
R = dvg.make_resolvent(n, arcs)
Us = [[(0,3),(0,4),(1,3),(1,4),(2,3),(2,4)], [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4)],
      [(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,4)], [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]]
B2 = PolynomialRing(QQ, ['z'] + ['a%d' % i for i in range(N)])
zz = B2.gens()[0]
for U in Us:
    U = frozenset(U)
    A = PolynomialRing(dvg.Kz, ['a%d' % i for i in range(N)], order='lex')
    al = A.gens(); memo = {}
    eqs = []; Ws = set()
    for (W, v) in P:
        eqs.append(dvg.mu_weighted((W - {v}) | F, U, al, A, memo) - A(R(W, v)) * dvg.mu_weighted(W | F, U, al, A, memo))
        Ws.add(W)
    I = A.ideal(eqs)
    for W in Ws:
        I = I.saturation(A.ideal([dvg.mu_weighted(W | F, U, al, A, memo)]))[0]
    G = I.groebner_basis()
    dim = I.dimension()
    report = []
    for var in range(N):
        J = I.elimination_ideal([al[j] for j in range(N) if j != var])
        g = J.gens()[0]
        # clear denominators and factor in Q[z, a_var]
        den = 1
        for c in g.coefficients():
            den = den.lcm(c.denominator()) if hasattr(den, 'lcm') else c.denominator()
        coeffs = g.dict()
        h = B2(0)
        for mon, c in coeffs.items():
            cc = c * den
            num = cc.numerator(); dd = cc.denominator()
            assert dd == 1 or dd.degree() == 0
            h += B2(num(zz)) / (dd if dd != 1 else 1) * B2.gens()[1 + var] ** mon[var]
        facs = h.factor()
        degs = [ (f.degree(B2.gens()[1 + var]), e) for f, e in facs if f.degree(B2.gens()[1 + var]) > 0]
        report.append((var, degs))
    print('U =', sorted(U), 'dim', dim, 'factor degrees of eliminants (var, [(deg,mult)]):', report)
