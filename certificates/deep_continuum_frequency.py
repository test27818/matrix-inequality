"""Integer/dyadic main transform versus the deep residual. See accompanying proof."""
import argparse
from fractions import Fraction as F
from math import factorial
from functools import lru_cache
from moderate_helmholtz import I,Q,ONE,ZERO
def ir(a,b):return I(I.q(a).lo,I.q(b).hi)
def series(z,n):
    # sinhc / coshc, <=9; through degree 82, tail < 10^-40.
    assert 0<=z.lo and z.hi<=9*Q
    if n==0:return (z.exp()+(-z).exp())/2
    w=z.sq(); term=I.q(F(1,factorial(n))); s=term
    for k in range(1,42):
        term=(term*w)/((2*k+n)*(2*k+n-1)); s+=term
    return s+I(0,I.q(F(1,10**40)).hi)
@lru_cache(maxsize=65536)
def pointcs(z):
    x=I.q(z);k=0
    while x.hi>Q: x=x/2; k+=1
    xx=x.sq(); co=si=ONE; a=b=ONE
    for j in range(1,15):
        a=-a*xx/((2*j-1)*(2*j)); b=-b*xx/((2*j)*(2*j+1))
        co+=a; si+=b
    co=I(co.lo-(a*xx/(29*30)).hi,co.hi);si=I(si.lo-(b*xx/(30*31)).hi,si.hi)
    for _ in range(k):
        si=si*co;co=2*co.sq()-1
    return co,si
def weights(N,w0,w1):
    w=(w0+w1)/2;delta=(w1-w0)/2
    sc=pointcs(w/F(2*N))[1]
    out=[]
    for i in range(N):
        pm=F(2*i+1,2*N)
        iv=1-pointcs(w*pm)[0]*sc
        err=I.q(pm*delta).hi
        out.append(I(max(0,iv.lo-err),min(2*Q,iv.hi+err)))
    return out
def alphas(c0,c1,N):
    c=ir(c0,c1);u=ir(0,F(1)/(24+c0))
    sc1=series(c,1);ch1=series(c,0);ac=series(c,2)
    A=[]
    for i in range(N):
        q=ir(1-F(i+1,N),1-F(i,N));z=q*c
        sh=series(z,1); ch=series(z,0)
        M2=q*(2*sh+ch);N2=2*ch+z.sq()*sh
        A.append((1+u)*M2-u*N2)
    high=sc1-u*ch1-((1+u)*(2*sc1+ch1)+u*(2*ch1+c.sq()*sc1))/32
    return A,ac,u,high
def main():
    pa=argparse.ArgumentParser();pa.add_argument('--mesh',type=int,default=64);pa.add_argument('--wgrid',type=int,default=4);pa.add_argument('--cgrid',type=int,default=8)
    aa=pa.parse_args();N=aa.mesh
    roots=[]
    for wi in range(aa.wgrid,32*aa.wgrid):
        w0=F(wi,aa.wgrid);w1=F(wi+1,aa.wgrid);roots.append((w0,w1,weights(N,w0,w1)))
    # Conservative ↑ roundings of --u 2 --s 4 bounds; use assert, not input!
    from deep_continuum_tables import MASS as mass, CURV as curvature
    # Alias tail, all derivatives at most 1/250, glued p>=1 with <1e-7.
    for k in range(9):
        worstlow=worstmid=worsthigh=10**8*Q
        bs=F(curvature[k],10**6);ms=F(mass[k],10**6)
        for j in range(aa.cgrid):
            c0=k+F(j,aa.cgrid);c1=c0+F(1,aa.cgrid)
            A,a,u,high=alphas(c0,c1,N)
            # omega<=1: near-origin scale; use rigorous integrated p² weights.
            s=ZERO
            for i in range(N):
                l=F(i,N);r=F(i+1,N);p2=(r**3-l**3)*N/3
                s+=A[i]*ir(p2*(F(1,2)-r*r/24),p2/2)/N
            low=s-a*u.sq()-I.q(ms+F(41,10000))
            worstlow=min(worstlow,low.lo)
            for w0,w1,ww in roots:
                # A and weights represented at common scale. Summation before outward division.
                val=sum((a0.lo*w.hi if a0.lo<0 else a0.lo*w.lo) for a0,w in zip(A,ww))//(N*Q)
                spike=a*(1+u)*I.q(w1*w1/((24+c0)**2+w1*w1))
                cost=min(2*bs,ms*w1*w1)+min(F(81,10000),F(41,10000)*w1*w1)
                v=val-spike.hi-I.q(cost).hi
                worstmid=min(worstmid,v)
            worsthigh=min(worsthigh,(high-I.q(2*bs+F(81,10000))).lo)
        print('c',k,k+1,'bound low mid high in millionths',*[x*10**6//Q for x in (worstlow,worstmid,worsthigh)],flush=True)
        assert min(worstlow,worstmid,worsthigh)>0
if __name__=='__main__': main()
