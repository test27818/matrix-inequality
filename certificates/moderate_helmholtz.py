#!/usr/bin/env python3
"""Independent exact outward-interval checker for moderate radial bridges.

Every endpoint is an integer multiple of 2**-112. No use of floating arithmetic
for certificate acceptance or parameter subdivision. See moderate-helmholtz.md
for the analytic interpretation. This is not the missing legacy C++ checker.
"""
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
from math import isqrt
import argparse

BITS = 112
Q = 1 << BITS


def ceildiv(a, b):
    return -((-a) // b)


@dataclass(frozen=True, slots=True)
class I:
    lo: int
    hi: int

    @staticmethod
    def q(x):
        if isinstance(x, I):
            return x
        x = Fraction(x)
        return I(x.numerator*Q//x.denominator, ceildiv(x.numerator*Q, x.denominator))

    def __add__(a, b):
        if isinstance(b, int):
            t = b*Q
            return I(a.lo+t, a.hi+t)
        b = I.q(b)
        return I(a.lo+b.lo, a.hi+b.hi)

    __radd__ = __add__

    def __neg__(a):
        return I(-a.hi, -a.lo)

    def __sub__(a, b):
        if isinstance(b, int):
            t = b*Q
            return I(a.lo-t, a.hi-t)
        b = I.q(b)
        return I(a.lo-b.hi, a.hi-b.lo)

    def __rsub__(a, b):
        return I.q(b)-a

    def __mul__(a, b):
        if isinstance(b, int):
            if b >= 0:
                return I(a.lo*b, a.hi*b)
            return I(a.hi*b, a.lo*b)
        b = I.q(b)
        if a.lo >= 0 and b.lo >= 0:
            return I(a.lo*b.lo//Q, ceildiv(a.hi*b.hi, Q))
        p = (a.lo*b.lo, a.lo*b.hi, a.hi*b.lo, a.hi*b.hi)
        return I(min(p)//Q, ceildiv(max(p), Q))

    __rmul__ = __mul__

    def __truediv__(a, b):
        if isinstance(b, int) and b > 0:
            return I(a.lo//b, ceildiv(a.hi, b))
        b = I.q(b)
        if b.lo <= 0 <= b.hi:
            raise ArithmeticError("division interval includes zero")
        if b.hi < 0:
            return -a/-b
        if a.lo >= 0:
            return I(a.lo*Q//b.hi, ceildiv(a.hi*Q, b.lo))
        if a.hi <= 0:
            return I(a.lo*Q//b.lo, ceildiv(a.hi*Q, b.hi))
        return I(a.lo*Q//b.lo, ceildiv(a.hi*Q, b.lo))

    def __rtruediv__(a, b):
        return I.q(b)/a

    def sq(a):
        lo = min(a.lo*a.lo, a.hi*a.hi) if not a.lo <= 0 <= a.hi else 0
        hi = max(a.lo*a.lo, a.hi*a.hi)
        return I(lo//Q, ceildiv(hi,Q))

    def sqrt(a):
        if a.lo < 0:
            raise ArithmeticError("sqrt interval has negative lower endpoint")
        lo, hi = isqrt(a.lo*Q), isqrt(a.hi*Q)
        return I(lo, hi + (hi*hi < a.hi*Q))

    def exp(a):
        l = _exp_at(a.lo)
        h = l if a.lo == a.hi else _exp_at(a.hi)
        return I(l.lo, h.hi)

    def point_mid(a):
        mid = (a.lo+a.hi)//2
        return I(mid,mid)

    def span(a):
        return a.hi-a.lo

    def split(a):
        mid = (a.lo+a.hi)//2
        if not a.lo < mid < a.hi:
            raise ArithmeticError("cannot split a single dyadic")
        return I(a.lo,mid),I(mid,a.hi)

ZERO, ONE = I(0,0), I(Q,Q)


@lru_cache(maxsize=8192)
def _exp_at(a):
    if a == 0:
        return ONE
    if a < 0:
        return ONE/_exp_at(-a)
    # Rescale into [0,1], where terms past 32 are bounded by a geometric series.
    k = max(0, a.bit_length()-BITS)
    x = I(a,a)/(1 << k)
    s = t = ONE
    for j in range(1, 33):
        t = (t*x)/j
        s = s+t
    tail = (t*x/33)/(1-x/34)
    w = I(s.lo,s.hi+tail.hi)
    for _ in range(k):
        w = w.sq()
    return w


@dataclass(frozen=True, slots=True)
class J:
    """Taylor lanes f, df, d²f/2, irrespective of box widths."""
    f: I
    d: I = ZERO
    dd: I = ZERO

    @staticmethod
    def q(x):
        return x if isinstance(x,J) else J(I.q(x))

    def __add__(a,b):
        b=J.q(b)
        return J(a.f+b.f,a.d+b.d,a.dd+b.dd)

    __radd__=__add__

    def __neg__(a):
        return J(-a.f,-a.d,-a.dd)

    def __sub__(a,b):
        return a+-J.q(b)

    def __rsub__(a,b):
        return J.q(b)+-a

    def __mul__(a,b):
        b=J.q(b)
        return J(a.f*b.f,a.f*b.d+a.d*b.f,a.f*b.dd+a.d*b.d+a.dd*b.f)

    __rmul__=__mul__

    def __truediv__(a,b):
        b=J.q(b)
        r=a.f/b.f
        s=(a.d-r*b.d)/b.f
        t=(a.dd-r*b.dd-s*b.d)/b.f
        return J(r,s,t)

    def __rtruediv__(a,b):
        return J.q(b)/a

    def e(a):
        r=a.f.exp()
        return J(r,r*a.d,r*(a.dd+a.d.sq()/2))

    def s(a):
        r=a.f.sqrt()
        u=2*r
        s=a.d/u
        return J(r,s,(a.dd-s.sq())/u)


def calc_h(t,b):
    """Ĥ=2*exp(t) ∫ exp(-v)/(1-exp(-v))²*sinh(κ v)/κ dv.

    q=1 alias and first 24 terms in b=κ². For 16≤t≤56, the q=1 b-tail
    is bounded by 2*e^56*(1/256)^24/(1-1/256) < 2^-100 using e<3.
    The other aliases multiply the integrand by ≤2*exp(-t)/(1-exp(-t))².
    """
    z = ONE
    moment = ONE
    cur = ONE
    accum = ZERO
    for j in range(1,48):
        z=z*t/j                 # t**j/j!
        moment=moment+z         # sum_{i=0}^j t**i/i!
        if j%2:
            accum=accum+cur*moment
            cur=cur*b
    leading=2*accum + I(0,Q >> 100)
    e=(-t).exp()
    denom=(1-e).sq()
    alias=2*e/denom
    return I(leading.lo,(leading*(1+alias)).hi),3*e/denom


def moment012(t,a):
    """Taylor magnitudes exp(t)*∫ v^i/i!*exp(-a v) dv, i=0,1,2."""
    inv=1/a
    e=((1-a)*t).exp()*inv
    t1=t+inv
    return (e,e*t1,e*(t.sq()/2+inv*t1))


def alias_hi(delta,moment):
    return I(0,(delta*moment).hi)


def tjet(t,p,data):
    """T̂_t=e^t(exp(-p t)/p-J_t), with the removable endpoint p=0."""
    h,delta=data
    a=moment012(t,2-p)
    d=moment012(t,1+p)
    A=tuple(alias_hi(delta,x) for x in a)
    B=tuple(alias_hi(delta,x) for x in d)
    return J((1-p)*h+a[0]-2*d[0]+A[0]-B[0],-h+a[1]+2*d[1]+A[1]+B[1],a[2]-2*d[2]+A[2]-B[2])


BETA = I.q(Fraction(29,32))


@dataclass(slots=True)
class Box:
    x: I
    y: I
    b: I

    def model(box):
        return Model(box.x,box.y,box.b)


class Model:
    __slots__=('x','y','c','m','xs','ms','ys','em','c2')
    def __init__(self,x,y,b):
        c=(y-x)/2
        m=(y+x)/2
        if not (16*Q <= x.lo <= x.hi <= 24*Q and 32*Q <= y.lo <= y.hi <= 56*Q
                and 0 <= b.lo <= b.hi <= Q//256):
            raise ValueError("radial outside supported domain")
        self.x,self.y,self.c,self.m=x,y,c,m
        self.xs,self.ms,self.ys=calc_h(x,b),calc_h(m,b),calc_h(y,b)
        self.em=(-c).exp()
        self.c2=c.sq()

    @staticmethod
    def W(t,p,T):
        return 1-p*((p-1)*J(t)).e()*T

    def low(self, pp):
        """Lower enclosure of e^x*(c²D−2max(c²D−D″,0))."""
        x,m,y,c=self.x,self.m,self.y,self.c
        p=J(pp,ONE)
        Tx,Tm,Ty=tjet(x,pp,self.xs),tjet(m,pp,self.ms),tjet(y,pp,self.ys)
        vx=self.W(x,p,Tx).s()
        vy=self.W(y,p,Ty).s()
        G=Tx/(1+vx)
        E=-J(self.em)*Tm+vx*((p-2)*J(c)).e()*Ty/(1+vy)
        decay=(-J(c)*p).e().f
        D=decay*G.f+E.f
        H=decay*(2*c*G.d-2*G.dd)+self.c2*E.f-2*E.dd
        L=self.c2*D
        return min(L.lo,(L-2*H).lo)

    def high(self):
        """Lower bound for composite tail plus the cutoff atom."""
        x,m,y,c=self.x,self.m,self.y,self.c
        p=J(BETA,ONE)
        Wx=self.W(x,p,tjet(x,BETA,self.xs))
        Wy=self.W(y,p,tjet(y,BETA,self.ys))
        # e^(x-beta*m)=e^((1-beta)*x-beta*c).
        ep=(J(x*(1-BETA))-p*J(c)-J(x)*(p-J(BETA))).e()
        V=(ep/p)*(Wx*Wy).s()
        return (2*V.d-self.c2*ep.f/(m*BETA)).lo


def mesh(power):
    if not 0<=power<=16:
        raise ValueError("mesh exponent out of range")
    K=1<<power
    width=BETA.hi//K
    return [I(j*width,(j+1)*width) for j in range(K)]


def one(box,nodes):
    """Outward lower sum; None only means that a denominator guard failed."""
    M=box.model()
    try:
        v=M.high()
        for p in nodes:
            v+=M.low(p)*p.span()//Q
        return v
    except ArithmeticError:
        return None


def exact_box(ranges):
    xi,xa,yi,ya,bi,ba = map(Fraction,ranges)
    if not (xi<=xa and yi<=ya and bi<=ba):
        raise ValueError("empty box")
    f=lambda i,a:I(I.q(i).lo,I.q(a).hi)
    return Box(f(xi,xa),f(yi,ya),f(bi,ba))


ROOTS = {
    'lower-left':(16,20,32,48,0,'1/256'),
    'lower-right':(20,24,32,48,0,'1/256'),
    'upper':(16,24,48,56,0,'1/256'),
}


def split_axis(box):
    x,y,b=box.x,box.y,box.b
    # Compare radii after rough kinematic and sinh(z)/z sensitivities.
    return max(range(3),key=lambda j: (x.span()*3,y.span()*2,b.span()*600)[j])


def children(box, axis):
    x,y,b=box.x,box.y,box.b
    if axis==0:
        return [Box(t,y,b) for t in x.split()]
    if axis==1:
        return [Box(x,t,b) for t in y.split()]
    return [Box(x,y,t) for t in b.split()]


def cover(box,nodes,maxdepth,limit,trace=None):
    todo=[(box,0)]
    low = None
    visits=accepted=0
    depth_max=0
    while todo:
        if visits>=limit:
            return False,visits,accepted,depth_max,low, 'limit',todo[-1][0]
        cell,d=todo.pop()
        visits+=1
        depth_max=max(depth_max,d)
        v=one(cell,nodes)
        if v is not None and v>0:
            low=v if low is None else min(low,v)
            accepted+=1
            if trace is not None:
                trace.append([0, v])
            continue
        if d>=maxdepth:
            return False,visits,accepted,depth_max,low,'depth',cell
        axis = split_axis(cell)
        L,R = children(cell,axis)
        if trace is not None:
            trace.append([1, axis, (L.x,L.y,L.b)[axis].hi])
        todo.append((R,d+1))
        todo.append((L,d+1))
    return True,visits,accepted,depth_max,low,'done',None


def simple(b):
    return None if b is None else tuple(str(Fraction(z,Q)) for i in (b.x,b.y,b.b) for z in (i.lo,i.hi))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',choices=list(ROOTS)+['all'])
    p.add_argument('--box',nargs=6,metavar=('X0','X1','Y0','Y1','B0','B1'))
    p.add_argument('--mesh',type=int,default=7,help='uniform quadrature mesh exponent')
    p.add_argument('--depth',type=int,default=24)
    p.add_argument('--limit',type=int,default=50000)
    a=p.parse_args()
    if bool(a.root)==bool(a.box):
        p.error('supply exactly one of --root and --box')
    roots = list(ROOTS) if a.root=='all' else [a.root] if a.root else ['custom']
    nodes=mesh(a.mesh)
    for name in roots:
        b=exact_box(a.box if a.box else ROOTS[name])
        ok,n,k,d,v,why,cell=cover(b,nodes,a.depth,a.limit)
        print(name, 'CERTIFIED' if ok else 'UNCERTIFIED',f'visits={n} leaves={k} depth={d}',
              f'lower={v}/{Q}' if v is not None else 'no accepted margin',why,simple(cell),flush=True)
        if not ok:
            return 1
    return 0

if __name__=='__main__':
    raise SystemExit(main())
