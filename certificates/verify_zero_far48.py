#!/usr/bin/env python3
"""Outward integer certificate for the true x∈[0,3], p∈[0,1/2] low-frequency jet.
No 1/x or 1/p; one fixed-point transcript at every continuous closed box.
"""
import argparse
from functools import lru_cache
from math import factorial
from fractions import Fraction as R
from verify_continuum_shallow import I,P,ival,frac,narrow,series

BA=I(0,P//256)
assert 2**26*P < factorial(25)

@lru_cache(None)
def exp_point(t):
    assert 0<=t<=2*P
    a=I(t)
    s=u=ival(1)
    for k in range(1,25):
        u=u*a/k
        s=s+u
    return I(s.l,s.h+1)

def ep(t):
    return I(exp_point(t.l).l,exp_point(t.h).h)

def bound(x):
    return max(-x.l,x.h)

class J:
    __slots__=("f","d","dd")
    def __init__(self,f,d=0,dd=0):
        self.f,self.d,self.dd=ival(f),ival(d),ival(dd)
    @staticmethod
    def of(x):
        return x if isinstance(x,J) else J(x)
    def __add__(self,z):
        z=J.of(z)
        return J(self.f+z.f,self.d+z.d,self.dd+z.dd)
    __radd__=__add__
    def __neg__(self):
        return J(-self.f,-self.d,-self.dd)
    def __sub__(self,z):
        return self+-J.of(z)
    def __rsub__(self,z):
        return J.of(z)+-self
    def __mul__(self,z):
        z=J.of(z)
        return J(self.f*z.f,self.f*z.d+self.d*z.f,self.f*z.dd+2*self.d*z.d+self.dd*z.f)
    __rmul__=__mul__
    def inv(self):
        r=1/self.f
        a=r.sq()
        return J(r,-self.d*a, (2*self.d.sq()*r-self.dd)*a)
    def __truediv__(self,z):
        return self*J.of(z).inv()
    def __rtruediv__(self,z):
        return J.of(z)*self.inv()
    def sqrt(self):
        r=self.f.sqrt()
        d=self.d/(2*r)
        return J(r,d,(self.dd-2*d.sq())/(2*r))

def ctail(q):
    """2 sum over j>=2; q jet is 1-p. Signs of the tail lanes are exact."""
    z=q*q
    n=48
    s=J(0)
    for j in range(2,n+1):
        s=s+2*j/((j*j-z)*(j*j-J(BA)))
    g=ival((n+1)**2)/((n+1)**2-BA)
    u=(n+1)**2-1
    s0=g/(n*(n+1))
    s1=s0*(frac(2,u))
    s2=s0*(frac(2,u)+frac(8,u*u))
    return s+J(I(0,s0.h),I(-s1.h,0),I(0,s2.h))

def profiles(pp,q,dv,n):
    prefix=[(ival(0),)*6]
    for j in range(n):
        v=I(j*dv.l,(j+1)*dv.h)
        z=v*q
        CV=series('C',v)
        fac=1/(1-BA*v.sq()/20)
        u0=v*series('E',z)/CV*fac
        u1=v.sq()*series('E1',z)/CV*fac
        u2=v.sq()*v*series('E2',z)/CV*fac
        h=1/ep(pp*v)
        pprev=prefix[-1]
        prefix.append(tuple(a+dv*b for a,b in zip(pprev,(u0,u1,u2,h,v*h,v.sq()*h))))
    return prefix

def run(nx,np):
    assert nx>0 and np>0 and P%nx==P%(2*np)==0
    dv=frac(1,nx); dphase=frac(1,2*np)
    maxim=[0,0,0]; where=[None,None,None]; fail=0; minW=P; ming0=99*P; slopes=[0,0]
    for k in range(np):
        pp=I(k*dphase.l,(k+1)*dphase.h)
        p=J(pp,1)
        q=1-p
        z=q*q
        Pq=q*(z-J(BA))
        regularpart=(2*z+2*q+1-J(BA))/((1+q)*(1-J(BA)))-Pq*ctail(q)
        pref=profiles(pp,q.f,dv,3*nx)
        for i in range(3*nx):
            x=I(i*dv.l,(i+1)*dv.h)
            cur=[I(pref[i][j].l,pref[i+1][j].h) for j in range(6)]
            # The p-derivative of the q-profile is negative.
            U=J(cur[0],-cur[1],cur[2])
            L=J(cur[3],-cur[4],cur[5])
            T=regularpart+Pq*U-L
            e=ep(pp*x)
            W=1-p*J(e,x*e,x.sq()*e)*T
            minW=min(minW,W.f.l)
            root=W.sqrt()
            G=T/(1+root)
            v=[bound(G.dd),bound(root.d),bound(root.dd)]
            for j in range(3):
                if v[j]>maxim[j]:
                    maxim[j]=v[j]
                    where[j]=(i,k)
            if k==0: ming0=min(ming0,G.f.l); slopes[0]=max(slopes[0],bound(G.d))
            if k==np-1: slopes[1]=max(slopes[1],bound(G.d))
            fail+=int(v[0]>4*P or v[1]>3*P or v[2]>21*P)
    print('rawMinW=',minW,'minG0=',ming0,'slopeBinMax=',slopes,'all divided by',P)
    fail+=int(minW<frac(9,64).h or ming0<P//16 or 5*slopes[0]>3*P or slopes[1]>P)
    print(f'radialBins={3*nx} spectralBins={np} boxes={3*nx*np} failures={fail}')
    print(f'maxAbs=[{maxim[0]}/{P},{maxim[1]}/{P},{maxim[2]}/{P}] cells={where}')
    return fail==0

def exp_neg_up(t):
    # Rational upper bound for exp(-t) by a positive partial sum of exp(t).
    den=a=R(1)
    for j in range(1,100):
        a*=R(t,j)
        den+=a
    return 1/den

def costs():
    y=R(48); w=(y+7)/2; cmin=R(45,2); cgmin=cmin/16
    E=(3*w*(1+y*y/8)*exp_neg_up(15*y/32)
       +R(3,2)*w*w*exp_neg_up(y/2)
       +3*(y+2)*(y*y/4+7*y/2+26)*exp_neg_up(11*y/16))/cgmin
    H=16*(2+4*R(17,15)+32/cmin)*exp_neg_up(12)
    A=2*exp_neg_up(cmin/2)/cgmin
    assert E<R(1,1000) and H<R(1,1000) and A<R(1,1000)
    residual=1-R(3,5)/cgmin-64/(cmin*cmin)-R(3,1000)
    assert residual>R(2,5)
    print('3 residual transform fees < 1/1000 each; normalized lower',residual,'> 2/5')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--radial-per-unit',type=int,default=32)
    ap.add_argument('--spectral-cells',type=int,default=32)
    a=ap.parse_args()
    costs()
    raise SystemExit(0 if run(a.radial_per_unit,a.spectral_cells) else 1)
