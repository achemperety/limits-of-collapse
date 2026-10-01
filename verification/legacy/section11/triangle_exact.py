"""Theorem 11.4, exactly: with t^3 - 3t = z, cycle weights alpha = z + t, and environment weights the roots of
(t^2 - 1) x^2 - 2 t x + 2, every frozen identity of the directed triangle holds, and the sextic is the eliminant."""
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField
R = PolynomialRing(QQ, ['t', 'x3', 'x4']); t, x3, x4 = R.gens()
z = t**3 - 3*t
alpha = z + t
# environment weights through their symmetric functions: s = x3 + x4, p = x3 x4
F = FractionField(PolynomialRing(QQ, 't')); T = F.gen()
zT = T**3 - 3*T; aT = zT + T
p = 2 / (T**2 - 1); s = T * p
m = {0: p, 1: aT * p - s, 2: aT**2 * p - 2 * aT * s + 2, 3: aT**3 * p - 3 * aT**2 * s + 6 * aT}
muP = {0: 1, 1: zT, 2: zT**2 - 1, 3: zT**3 - 2 * zT}
print('m_k = p * mu(P_k)(z) for k = 0..3:', all(m[k] == p * muP[k] for k in range(4)))
# frozen identities: c(interval of length j) = mu(P_{3-j}) / mu(P_3)
print('frozen identities:', all(m[3 - j] / m[3] == F(muP[3 - j]) / muP[3] for j in range(4)))
# elimination: resultant in t of (t^2-1)a^2 - 2 t a + 2 and t^3 - 3t - z gives the sextic
P2 = PolynomialRing(QQ, ['t', 'a', 'Z']); tt, a, Z = P2.gens()
res = ((tt**2 - 1) * a**2 - 2 * tt * a + 2).resultant(tt**3 - 3 * tt - Z, tt)
sextic = (Z**2 - 4) * a**6 + 12 * a**4 + 4 * Z * a**3 - 12 * a**2 + 8
print('resultant is a constant multiple of the sextic:', (res * sextic.lc() - sextic * res.lc()) == 0, '| resultant =', res.factor())
print('sextic = (z a^3 + 2)^2 + 4 (1 - a^2)^3:', sextic == (Z * a**3 + 2)**2 + 4 * (1 - a**2)**3)
