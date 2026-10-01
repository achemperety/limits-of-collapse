import sympy as sp
t,a,zz=sp.symbols('t a z')
z=t**3-3*t
al=z+t
p=2/(t**2-1); s=2*t/(t**2-1)
m0=p; m1=al*p-s; m2=al**2*p-2*al*s+2; m3=al**3*p-3*al**2*s+6*al
print('m0/m1 - 1/z:', sp.simplify(m0/m1-1/z))
print('m1/m2 - z/(z^2-1):', sp.simplify(m1/m2-z/(z**2-1)))
print('m2/m3 - (z^2-1)/(z^3-2z):', sp.simplify(m2/m3-(z**2-1)/(z**3-2*z)))
sext=(zz**2-4)*a**6+12*a**4+4*zz*a**3-12*a**2+8
print('sextic identity:', sp.expand(sext-((zz*a**3+2)**2+4*(1-a**2)**3)))
print('sextic irreducible over Q(z):', sp.Poly(sext,a,zz).is_irreducible)
# environment weight x satisfies the sextic with z=t^3-3t ?
x=sp.symbols('x')
quad=(t**2-1)*x**2-2*t*x+2
res=sp.resultant(quad, sext.subs(zz,z).subs(a,x), x)
print('resultant(quad, sextic) == 0:', sp.simplify(res)==0)
# cubic t^3-3t-z irreducible over Q(z)
print('cubic irreducible:', sp.Poly(t**3-3*t-zz,t,zz).is_irreducible)
