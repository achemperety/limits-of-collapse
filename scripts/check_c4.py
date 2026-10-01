import sympy as sp, itertools
z,p,q=sp.symbols('z p q', real=True)
y=sp.symbols('y0:4', real=True)
def mu_multi(n, edges, w):
    # multivariate matching polynomial sum_M (-1)^|M| prod_{unmatched} w
    tot=0
    E=list(edges)
    for r in range(0,3):
        for M in itertools.combinations(E,r):
            vs=[v for e in M for v in e]
            if len(set(vs))<len(vs): continue
            term=(-1)**r
            for v in range(n):
                if v not in vs: term*=w[v]
            tot+=term
    return sp.expand(tot)
Q=mu_multi(4,[(0,1),(1,2),(2,3),(3,0)],y)
muP=lambda k: [1,z,z**2-1,z**3-2*z,z**4-3*z**2+1][k]
P=sp.expand(Q+(z**2-1)+(p-1)*(y[0]-z)*(y[2]-z)+(q-1)*(y[1]-z)*(y[3]-z))
# check Taylor coefficients of P at y=z: for intervals = mu(P_{4-|T|}), for T=empty = mu(P_4)
sub={yy:z for yy in y}
def coeff(T):
    e=P
    for i in T: e=sp.diff(e,y[i])
    return sp.expand(e.subs(sub))
ints=[(),(0,),(1,),(2,),(3,),(0,1),(1,2),(2,3),(3,0),(0,1,2),(1,2,3),(2,3,0),(3,0,1),(0,1,2,3)]
print('interval coeffs ok:', all(sp.simplify(coeff(T)-muP(4-len(T)))==0 for T in ints))
print('coeff {0,2}:', coeff((0,2)), ' coeff {1,3}:', coeff((1,3)))
def Delta(P,i,j):
    return sp.expand(sp.diff(P,y[i])*sp.diff(P,y[j])-P*sp.diff(P,y[i],y[j]))
d02=sp.factor(Delta(P,0,2).subs({y[1]:z,y[3]:z}))
print('Delta02 at y1=y3=z:', d02, '  target 1-p*mu(P4):', sp.simplify(d02-(1-p*muP(4))))
d13=sp.factor(Delta(P,1,3).subs({y[0]:0,y[2]:0}))
print('Delta13 at y0=y2=0:', d13, ' diff:', sp.simplify(d13-(1-q)*(z**2*p+1)))
d12=sp.expand(Delta(P,1,2).subs({y[0]:-1,y[3]:1}))
print('Delta12 at y0=-1,y3=1:', d12)
print('  diff from claim:', sp.simplify(d12-(4-z**2+z**2*(p+q)+z*(p-q)+p*q*(z**2-1))))
# bound: maximize over box p,q in [-1/z^2, 1/mu(P4)] at z>=3: bilinear -> vertices
lo=-1/z**2; hi=1/muP(4)
for zz in [3,4,10,100]:
    vals=[ (d12).subs({p:a,q:b,z:zz}) for a in (lo,hi) for b in (lo,hi)]
    vals=[sp.nsimplify(v.subs(z,zz)) for v in vals]
    print(zz, max(float(v) for v in vals), 8-zz**2)
