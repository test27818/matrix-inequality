"""Exact dyadic enclosures of the continuum boundary-layer square and first B alias.
No floating arithmetic in acceptance. Analytic derivation: deep-continuum-closure.md.
"""
from fractions import Fraction as F
from math import factorial
from functools import lru_cache
from moderate_helmholtz import I,J,Q,ZERO,ONE
def ir(a,b): return I(I.q(a).lo,I.q(b).hi)
def absmax(x):return max(-x.lo,x.hi)
def sym(a):return I(-a,a)
def rat(x):return F(x,Q)
def vabs(j):return [F(absmax(v),Q) for v in (j.f,j.d,j.dd)]
def jerr(j,eps): return J(j.f+eps,j.d+eps,j.dd+eps)
def ec(t,n):
    # E_n(t)=sum t^{2j}/(2j+n)!  (n=2,3), t<=2.6 and |dt/dT|<=1.3.
    # Coefficient tail past index 16 in f,f',f''/2 <= 10^-22.
    assert 0<=t.f.lo and t.f.hi<=3*Q and absmax(t.d)<=Q*13//10 and t.dd==ZERO
    w=t*t
    z=J.q(0)
    for k in range(16,-1,-1):
        z=z*w+F(1,factorial(2*k+n))
    ep=I(0,I.q(F(1,10**22)).hi)
    return jerr(z,ep)
@lru_cache(maxsize=16384)
def ke(u,r,Tf):
    T=J(Tf,ONE); qf=u*Tf
    assert qf.lo<=Q
    q=J(I(qf.lo,min(qf.hi,Q)),u)
    p=1-q
    t=T*r
    if Tf.hi<=2*Q:
        a=1/r*u; E=ec(t,3); C=ec(t,2); ecc=(-t).e()
        # W/t^3. All nonleading aliases have negligible normalized T-jets.
        W=jerr(2*ecc*(E+C*a+a.sq()+a.sq()*a)/(1+q),sym(I.q(F(1,10**7)).hi))
        numer=(C/r+u/r.sq()+u.sq()/(r.sq()*r)-T*E)/(1+q)
        numer=jerr(numer,sym(I.q(F(1,10**6)).hi))
        return (T/2-t).e()*(r*r.sqrt())*numer/W.s()
    e=(-t).e()
    W=jerr(1-p*(e*e)/(1+q)-2*p*(t+q)*e,sym(I.q(F(1,10**7)).hi))
    if W.f.lo<=0: raise ArithmeticError(('wbad',rat(Tf.lo),rat(Tf.hi),rat(r.lo),rat(r.hi)))
    lead=(T/2-t).e()*T*(e/(1+q)+t-p)
    # b0=m exp(-v+T/2). Bounds v>=24 and v>=m-9; m>=T.
    # In each T cell use minimum of two upper bounds, safely for all c<=9.
    M0=max(24*Q,Tf.lo); M0=I(M0,M0)
    bd0=(33*(Tf/2-24).exp()).hi
    bd1=(M0*(Tf/2+9-M0).exp()).hi
    bd=I(0,min(bd0,bd1))
    coeff=J(bd,bd/2,bd/8)
    j1=coeff*(-1/(1+p)+(1+2*q)/(2+q)*e*e+((1+q)*t-p*q/2)*e)
    small=sym((bd*100*I.q(-24).exp()).hi)
    return jerr(lead+j1,small)/W.s()
def layer(c0,c1,u0,u1,t0,t1,ns):
    u=ir(u0,u1);Tf=ir(t0,t1)
    if (u*Tf).lo>Q:return (F(0),F(0))
    u=I(u.lo,min(u.hi,(1/I(Tf.lo,Tf.lo)).hi)) if Tf.lo else u
    hm=I((u*I.q(c0)).lo,min((u*I.q(c1)).hi, I.q(F(c1,24+c1)).hi))
    T=J(Tf,ONE)
    qf=u*Tf;q=J(I(qf.lo,min(qf.hi,Q)),u)
    accum=J.q(0);Bs=J.q(0)
    for k in range(ns):
        sr=ir(-1+F(2*k,ns),-1+F(2*k+2,ns))
        d=hm*sr; r=1+d
        accum+=ke(u,r,Tf)
        Bs+=T*(T*(d-F(1,2))).e()
    av=accum/ns
    B=Bs/ns
    if Tf.hi<=2*Q:
        Z=(1-q)*T*T*T*av*av
    else:
        Z=(1-q)*av*av
    Bj=B*B/(1+q)
    return F(absmax(Z.f)+absmax(Bj.f),Q),F(2*(absmax(Z.dd)+absmax(Bj.dd)),Q)
def bins(a,b,num):
    for i in range(num): yield a+(b-a)*i/num,a+(b-a)*(i+1)/num
def box(c0,c1,uz=4,ns=4,scale=1):
    umax=F(1,24+c0)
    Tpart=list(bins(F(0),F(2),32*scale))+list(bins(F(2),F(12),80*scale))+list(bins(F(12),F(36),96*scale))+list(bins(F(36),F(96),120*scale))
    vals=[]
    for u0,u1 in bins(F(0),umax,uz):
        a=b=F(0)
        for t0,t1 in Tpart:
            A,B=layer(c0,c1,u0,u1,t0,t1,ns)
            a+=A*(t1-t0);b+=B*(t1-t0)
        vals.append((u1**4*(a+F(1,10**8)),u1**2*(b+F(1,10**8))))
    return tuple(max(a[j] for a in vals) for j in (0,1))
if __name__=='__main__':
    import argparse
    par=argparse.ArgumentParser();par.add_argument('left',type=int,nargs='?',default=0);par.add_argument('right',type=int,nargs='?',default=9);par.add_argument('--u',type=int,default=2);par.add_argument('--s',type=int,default=4);par.add_argument('--scale',type=int,default=1)
    args=par.parse_args()
    assert 0<=args.left<args.right<=9
    for k in range(args.left,args.right):
        A,B=box(k,k+1,args.u,args.s,args.scale)
        print('c',k,k+1,'mass <=',-(-A*10**6//1),'/1000000','integrated second derivative <=',-(-B*10**6//1),'/1000000',flush=True)
        from deep_continuum_tables import MASS,CURV
        if A*10**6>MASS[k] or B*10**6>CURV[k]: raise ArithmeticError('candidate table bound failed')
