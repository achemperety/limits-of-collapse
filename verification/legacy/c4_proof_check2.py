from sage.all__sagemath_singular import PolynomialRing, QQ
R = PolynomialRing(QQ, ['z', 'p', 'q', 'y0', 'y1', 'y2', 'y3'])
z, p, q, y0, y1, y2, y3 = R.gens()
y = [y0, y1, y2, y3]
u = p - 1; v = q - 1
Q = y0*y1*y2*y3 - y0*y1 - y1*y2 - y2*y3 - y3*y0 + 2
P = Q + (z**2 - 1) + u*(y0 - z)*(y2 - z) + v*(y1 - z)*(y3 - z)
def Delta(i, j):
    return P.derivative(y[i]) * P.derivative(y[j]) - P * P.derivative(y[i]).derivative(y[j])
print('Delta_12(y0=-1,y3=1) =', Delta(1, 2).subs({y0: -1, y3: 1}))
print('Delta_13(y0=y2=0)    =', Delta(1, 3).subs({y0: 0, y2: 0}).factor())
print('Delta_02(y1=y3=z)    =', Delta(0, 2).subs({y1: z, y3: z}))
muP4 = z**4 - 3*z**2 + 1
# Rayleigh at y = z for the non-adjacent pair {0,2}:  c0 c2 - c_empty c02 >= 0 in Taylor coefficients
T0 = P.derivative(y0).subs({y0: z, y1: z, y2: z, y3: z}); T2 = P.derivative(y2).subs({y0: z, y1: z, y2: z, y3: z})
T02 = P.derivative(y0).derivative(y2).subs({y0: z, y1: z, y2: z, y3: z}); Te = P.subs({y0: z, y1: z, y2: z, y3: z})
print('Rayleigh {0,2} at y=z : T0*T2 - Te*T02 =', (T0 * T2 - Te * T02))
