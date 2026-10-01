"""Mean-field environments (every cycle vertex has the same neighbours N in F) with U[C] empty (kappa = 0) or
complete (kappa = 1). The frozen identities say that e^{-alpha x + kappa x^2/2} * lambda * sum_k mu(P_k) x^k/k!
agrees with a polynomial of degree |N| up to x^m. With z = w + 1/w the coefficient of x^r/r! is
G_r(alpha) = lambda [w H_r(w - alpha) - w^{-1} H_r(w^{-1} - alpha)] / (w - w^{-1}),  H_r(A) = r! [x^r] e^{A x + kappa x^2/2}.
|N| <= m - 2 needs G_{m-1}(alpha) = G_m(alpha) = 0. We test whether the resultant in alpha vanishes identically in w."""
import sympy as sp
w, al, x = sp.symbols('w alpha x')
for kappa in (0, 1):
    for m in range(3, 11):
        def G(r):
            A1, A2 = w - al, 1 / w - al
            H = lambda A: sp.factorial(r) * sp.series(sp.exp(A * x + kappa * x ** 2 / 2), x, 0, r + 1).removeO().coeff(x, r)
            return sp.together(sp.expand(w * H(A1) - H(A2) / w))
        g1, g2 = sp.numer(G(m - 1)), sp.numer(G(m))
        res = sp.resultant(sp.Poly(g1, al), sp.Poly(g2, al))
        res = sp.factor(sp.simplify(res))
        print('kappa=%d m=%2d  resultant vanishes identically: %s   %s' % (kappa, m, res == 0, str(res)[:90]), flush=True)
