"""Theorem 11.18: the Jensen polynomial g_3 of the mean-field environment, and its two non-real roots.
Symbolic identities, then a numerical illustration at real z = 3 on the complex family of Theorem 11.5."""
import sympy as sp
w, al, t, lam = sp.symbols('w alpha t lambda')
A, B = w - al, 1 / w - al
F = lambda k: (w * A ** k - B ** k / w) / (w - 1 / w)
g3 = sum(sp.binomial(3, k) * lam * F(k) * t ** k for k in range(4))
closed = lam * (w * (1 + A * t) ** 3 - (1 + B * t) ** 3 / w) / (w - 1 / w)
print('g3 closed form:', sp.simplify(g3 - closed) == 0)
zeta = sp.Symbol('zeta')
tz = (zeta - 1) / (A - zeta * B)
val = sp.simplify((w * (1 + A * tz) ** 3 - (1 + B * tz) ** 3 / w) * (A - zeta * B) ** 3 / (A - B) ** 3)
print('root identity  w(1+At)^3 - w^-1(1+Bt)^3 = (A-B)^3 (w zeta^3 - 1/w)/(A - zeta B)^3:', sp.simplify(val - (w * zeta ** 3 - 1 / w)) == 0)
# numerical illustration: real z = 3, any real alpha -> g3 has two non-real roots
import random
wv = (3 + sp.sqrt(5)) / 2
for _ in range(5):
    a = sp.Rational(random.randint(-500, 500), 37)
    roots = sp.Poly(sp.expand(closed.subs({w: wv, al: a, lam: 1}) * (wv - 1 / wv)), t).nroots()
    print('  z = 3, alpha = %s: roots of g3 =' % a, [complex(r) for r in roots])
