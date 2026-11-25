"""Exact outward-rounded interval certificate for the continuum midpoint defect.
 The parameter q here is 1 - frequency index p, NOT the finite-order kappa.
 Usage: python verify_continuum_shallow.py [edge_numerator=2] [bins_per_unit=16] [phase_bins=32] [edge_denominator=1] [levels=4] [first_depth_cap_numerator] [first_depth_cap_denominator=1]
 Only integer operations / isqrt in bounds; timings are not part of the proof.
 See ../direct-scalar-route.md for the identity verified and mathematical derivation.
"""
from math import factorial, isqrt
from fractions import Fraction
from functools import lru_cache
import sys
P=1<<56
def ceildiv(a,b):return -((-a)//b)
class I:
    __slots__=('l','h')
    def __init__(self,l,h=None): self.l=l;self.h=l if h is None else h;assert self.l<=self.h
    def __add__(x,y):
        y=ival(y);return I(x.l+y.l,x.h+y.h)
    __radd__=__add__
    def __neg__(x):return I(-x.h,-x.l)
    def __sub__(x,y):return x+-ival(y)
    def __rsub__(x,y):return ival(y)+-x
    def __mul__(x,y):
        y=ival(y);a,b,c,d=x.l*y.l,x.l*y.h,x.h*y.l,x.h*y.h
        return I(min(a,b,c,d)//P,ceildiv(max(a,b,c,d),P))
    __rmul__=__mul__
    def __truediv__(x,y):
        y=ival(y);assert y.l>0
        return x*I(P*P//y.h,ceildiv(P*P,y.l))
    def __rtruediv__(x,y):return ival(y)/x
    def sqrt(x):assert x.l>=0;return I(isqrt(P*x.l),isqrt(P*x.h)+1)
    def sq(x):
        if x.h<0:return (-x).sq()
        if x.l<0:return I(0,ceildiv(max(-x.l,x.h)**2,P))
        return x*x
    def __repr__(x):return f'[{x.l}/{P},{x.h}/{P}]'
def ival(x):return x if isinstance(x,I) else I(P*x)
def frac(a,b):return I(P*a//b,ceildiv(P*a,b))
def narrow(x,y):
    return I(max(x.l,y.l),min(x.h,y.h))
def cp(f,n):
    # Exact coefficient = u/v! at degree n.
    if f=='cosh' and n%2==0:return 1,n
    if f=='sinhc' and n%2==0:return 1,n+1
    if f=='C' and n%2==0:return 1,n+2
    if f=='S' and n%2==0:return n,n+2
    if f=='E' and n%2==0:return 1,n+3
    if f=='E1' and n%2==1:return n+1,n+4
    if f=='E2' and n%2==0:return (n+2)*(n+1),n+5
    return 0,0
assert 2*4**49*P<factorial(49) # Taylor tail of each function at all 0<=v<=4
@lru_cache(None)
def series_at(f,x):
    assert 0<=x<=4*P
    l=h=0
    for n in range(49):
        u,k=cp(f,n)
        if not u:continue
        num=u*x**n*P;den=factorial(k)*P**n
        l+=num//den;h+=ceildiv(num,den)
    return l,h+1
def series(f,x):return I(series_at(f,x.l)[0],series_at(f,x.h)[1])
def phi(q):
    # phi(z)=(1-z)V, z=q^2. All tail bounds beyond j=n follow from
    # sum 1/[j(j^2-1)]=1/(2 n(n+1)) for j>n.
    z=q*q;n=48
    f=ival(2);fz=ival(0);fzz=ival(0)
    for j in range(2,n+1):
        d=j*j-z
        f=f+(2*(1-z))/(j*d)
        fz=fz-ival(2*(j*j-1))/(j*d*d)
        fzz=fzz-ival(4*(j*j-1))/(j*d*d*d)
    f=f+I(0,frac(1,n*(n+1)).h)
    fz=fz+I(-frac(1,n*(n+1)).h,0)
    fzz=fzz+I(-frac(2,n*(n+1)*((n+1)**2-1)).h,0)
    return (f,2*q*fz,2*fz+4*z*fzz)
def fvals(q,v):
    t=q*v
    C=series('C',v);S=series('S',v)
    H=(series('cosh',t)-(S/C)*series('sinhc',t))/C
    return v*series('E',t)/C,(v*v)*series('E1',t)/C,(v*v*v)*series('E2',t)/C,H
def rt(f):
    s,s1,s2=f; r=s.sqrt();return r,s1/(2*r),s2/(2*r)-s1.sq()/(4*r*r*r)
def run(num,bpu,nbq,den=1, targets=None,cap_num=None,cap_den=1):
    assert 0<num<=4*den and num*bpu%den==0 and P%bpu==0 and P%nbq==0
    n=num*bpu//den;dv=frac(1,bpu);dq=frac(1,nbq)
    best=10**100;which=None;failed=[]; passed_min=10**100; visited=0
    active=None if targets is None else {v[2] for v in targets}
    maxi=n if cap_num is None else min(n,cap_num*bpu//cap_den)
    assert cap_num is None or cap_num*bpu%cap_den==0
    byk={}
    if targets is not None:
        for i,m,k in targets:byk.setdefault(k,[]).append((i,m))
    for k in range(nbq):
        if active is not None and k not in active: continue
        q=I(k*dq.l,(k+1)*dq.h); z=q*q
        ph0,ph1,ph2=phi(q)
        f=[fvals(q,I(i*dv.l,(i+1)*dv.h)) for i in range(n)]
        su=[(ival(0),ival(0),ival(0))]
        for row in f:su.append(tuple(su[-1][j]+dv*row[j] for j in range(3)))
        roots=[]
        for i in range(n):
            U,V,X=(I(su[i][j].l,su[i+1][j].h) for j in range(3))
            roots.append(rt((ph0-(1-z)*U,ph1+2*q*U-(1-z)*V,ph2+2*U+4*q*V-(1-z)*X)))
        tabs=[f]
        while 2**len(tabs)<=n:
            last=tabs[-1];jump=2**(len(tabs)-1)
            tabs.append([tuple(I(min(last[i][j].l,last[i+jump][j].l),max(last[i][j].h,last[i+jump][j].h)) for j in range(4)) for i in range(len(last)-jump)])
        def rmq(lo,hi):
            j=(hi-lo+1).bit_length()-1;left=tabs[j][lo];right=tabs[j][hi+1-(1<<j)]
            return [I(min(x.l,y.l),max(x.h,y.h)) for x,y in zip(left,right)]
        pre=[[0,0,0,0,0,0]]
        for row in f:pre.append([pre[-1][j]+(row[j//2].l if j%2==0 else row[j//2].h) for j in range(6)])
        pairs=((i,m) for i in range(maxi) for m in range(i,n)) if targets is None else sorted(byk[k])
        for i,m in pairs:
                fs=rmq(i,m)
                visited+=1
                P0=fs[3]/8
                W=fs[:3]
                if m>=i+2:
                    for j in range(3):
                        # All intervening full radial cells, plus subsets of the two end cells
                        l=pre[m][2*j]-pre[i+1][2*j];h=pre[m][2*j+1]-pre[i+1][2*j+1]
                        W[j]=narrow(W[j],I(l//(m-i+1),ceildiv(h+f[i][j].h+f[m][j].h,m-i-1)))
                # F(q,v)<=v/3 implies W<=(A+B)/6
                W[0]=narrow(W[0],I(0,frac(i+m+2,6*bpu).h))
                S,S1,S2=(roots[i][j]+roots[m][j] for j in range(3))
                L=1/(2*S*S); l1=-2*S1/S; l2=6*(S1/S).sq()-2*S2/S
                def lower(p0,weights):
                    w,v,x=weights
                    Y=L*w.sq(); Y1=L*(2*w*v+w.sq()*l1);Y2=L*(2*v.sq()+2*w*x+4*w*v*l1+w.sq()*l2)
                    return p0+(6-20*z)*Y+(6*q-10*q*z)*Y1+(z-z*z)*Y2
                val=lower(P0,W)
                if val.l<=0 and m>i:
                    # W is an ordinary uniform average under v=A+(B-A)y.
                    # P0 is the average of H against min(y,1-y)/2.
                    M=16; wt=[ival(0),ival(0),ival(0)];pp=ival(0)
                    for slot in range(M):
                        lo=i+(m-i)*slot//M
                        hi=i+ceildiv((m-i)*(slot+1),M)
                        box=rmq(lo,hi)
                        for j in range(3):wt[j]=wt[j]+box[j]/M
                        hh=box[3]
                        pp=pp+hh*frac(2*min(slot,M-1-slot)+1,4*M*M)
                    val=lower(narrow(P0,pp),[narrow(x,y) for x,y in zip(W,wt)])
                if val.l<best: best,which=val.l,(i,m,k)
                if val.l<=0:failed.append((i,m,k))
                else:passed_min=min(val.l,passed_min)
    print('edge='+str(Fraction(num,den)),'radial_bins='+str(n),'phase_bins='+str(nbq),
          'visited='+str(visited),'unresolved='+str(len(failed)),
          'certified_lower='+str(best)+'/'+str(P),'worst_cell='+str(which),flush=True)
    if failed:
        print('certified_strips_in_original_domain first_depth<='+str(Fraction(min(w[0] for w in failed),bpu))+' OR second_depth<='+str(Fraction(min(w[1] for w in failed),bpu)),flush=True)
    return failed,passed_min
def check(num=2,bpu=16,nbq=32,den=1,levels=4,cap_num=None,cap_den=1):
    if cap_num is not None:print('first radial depth capped at',str(Fraction(cap_num,cap_den)),flush=True)
    targets=None; minbound=10**100
    for _ in range(levels):
        f,p=run(num,bpu,nbq,den,targets,cap_num,cap_den);minbound=min(minbound,p)
        if not f:
            print('PASS: all original boxes certified; minimum across accepted boxes',str(minbound)+'/'+str(P),flush=True);return True
        targets={(2*i+x,2*m+y,2*k+z) for i,m,k in f for x in [0,1] for y in [0,1] for z in [0,1] if 2*i+x<=2*m+y}
        bpu*=2;nbq*=2
    print('INCOMPLETE; no conclusion on unresolved boxes');return False
if __name__=='__main__':
    args=list(map(int,sys.argv[1:]));sys.exit(0 if check(*args) else 1)
