#!/usr/bin/env python3
"""Independent exact reconstruction of [3,16] × [32,48] Helmholtz windows.

Reuses the checked integer interval/jet arithmetic, but retains up to three
aliases for shallow depths and places the high-band cutoff at 3/4.
The high-band formula is exact, not a sampled-frequency assertion.
"""
from fractions import Fraction
import argparse
import json
from pathlib import Path
from moderate_helmholtz import Q, I, J, ZERO, ONE, exact_box, moment012, alias_hi, calc_h, tjet
from moderate_helmholtz import Box, children, split_axis, simple

BETA = I.q(Fraction(3,4))


def hyperbolic_column(t,b):
    """Return fixed positive alias enclosure, tail factor, and retained count."""
    if t.lo >= 10*Q:
        return (*calc_h(t,b),1)
    # Retain q=1,2,3; for q>=1 and j>=24, the b-series tail in a
    # single *weighted* q alias is at most 2 e^56 b^24/(1-b)<2^-100.
    K = 3
    bands = []
    for q in range(1,K+1):
        term=bpow=ONE
        moment=ONE/q
        s=ZERO
        for n in range(1,48):
            term=term*t/n
            moment=(moment+term)/q
            if n%2:
                s=s+bpow*moment
                bpow=bpow*b
        tail=I(0,Q>>100)
        bands.append(2*q*((1-q)*t).exp()*s+tail)
    ee=(-t).exp()
    denom=(1-ee).sq()
    next_decay=(-K*t).exp()
    u=alias_hi((K+1)*next_decay/denom,bands[0])
    return (sum(bands,ZERO)+u,(K+2)*next_decay/denom,K)


def expanded_tjet(t, p, data):
    h,δ,K = data
    if K == 1:
        return tjet(t,p,(h,δ))
    a = [ZERO]*3
    d = [ZERO]*3
    for j in range(1,K+1):
        Ai=moment012(t,j+1-p)
        Bi=moment012(t,j+p)
        if j==1:
            A=tuple(alias_hi(δ,v) for v in Ai)
            B=tuple(alias_hi(δ,v) for v in Bi)
        for k in range(3):
            a[k]=a[k]+j*Ai[k]
            d[k]=d[k]+(j+1)*Bi[k]
    return J((1-p)*h+a[0]-d[0]+A[0]-B[0],-h+a[1]+d[1]+A[1]+B[1],a[2]-d[2]+A[2]-B[2])


class Model:
    __slots__=('x','y','c','m','xs','ms','ys','em','c2')
    def __init__(self,x,y,b,single_alias=False):
        c=(y-x)/2
        m=(y+x)/2
        if not (3*Q <= x.lo <= x.hi <= 16*Q and 32*Q <= y.lo <= y.hi <= 48*Q
                and 0 <= b.lo <= b.hi <= Q//256):
            raise ValueError('outside reconstructed root shapes')
        self.x,self.y,self.c,self.m=x,y,c,m
        hcol=(lambda t: (*calc_h(t,b),1)) if single_alias else (lambda t: hyperbolic_column(t,b))
        self.xs,self.ms,self.ys=map(hcol,(x,m,y))
        self.em=(-c).exp()
        self.c2=c.sq()

    @staticmethod
    def W(t,p,T):
        return 1-p*((p-1)*J(t)).e()*T

    def low(self,pp):
        x,m,y,c=self.x,self.m,self.y,self.c
        p=J(pp,ONE)
        tx=expanded_tjet(x,pp,self.xs)
        tm=expanded_tjet(m,pp,self.ms)
        ty=expanded_tjet(y,pp,self.ys)
        sx=self.W(x,p,tx).s()
        sy=self.W(y,p,ty).s()
        G=tx/(1+sx)
        E=-J(self.em)*tm+sx*((p-2)*J(c)).e()*ty/(1+sy)
        decay=(-J(c)*p).e().f
        D=decay*G.f+E.f
        H=decay*(2*c*G.d-2*G.dd)+self.c2*E.f-2*E.dd
        L=self.c2*D
        return min(L.lo,(L-2*H).lo)

    def high(self):
        x,m,y,c=self.x,self.m,self.y,self.c
        p=J(BETA,ONE)
        Wx=self.W(x,p,expanded_tjet(x,BETA,self.xs))
        Wy=self.W(y,p,expanded_tjet(y,BETA,self.ys))
        ep=(J(x*(1-BETA))-p*J(c)-J(x)*(p-J(BETA))).e()
        V=(ep/p)*(Wx*Wy).s()
        return (2*V.d-self.c2*ep.f/(m*BETA)).lo


def one(box,nodes,single_alias=False):
    M=Model(box.x,box.y,box.b,single_alias)
    try:
        return M.high()+sum(M.low(p)*p.span()//Q for p in nodes)
    except ArithmeticError:
        return None


def cover(box,power,maxdepth,limit,trace=None,single_alias=False):
    width=BETA.hi//(1<<power)
    nodes=[I(j*width,(j+1)*width) for j in range(1<<power)]
    todo=[(box,0)]
    lo=None
    visits=accepted=depth=0
    while todo:
        if visits>=limit:
            return False,visits,accepted,depth,lo,'limit',todo[-1][0]
        cell,d=todo.pop()
        visits+=1
        depth=max(depth,d)
        v=one(cell,nodes,single_alias)
        if v is not None and v>0:
            accepted+=1
            lo=v if lo is None else min(lo,v)
            if trace is not None:
                trace.append([0,v])
            continue
        if d>=maxdepth:
            return False,visits,accepted,depth,lo,'depth',cell
        axis=split_axis(cell)
        L,R=children(cell,axis)
        if trace is not None:
            trace.append([1,axis,(L.x,L.y,L.b)[axis].hi])
        todo.append((R,d+1))
        todo.append((L,d+1))
    return True,visits,accepted,depth,lo,'done',None


ROOTS={name:(x0,x1,32,48,0,'1/256') for name,x0,x1 in [('left',3,9),('mid',9,12),('right',12,16)]}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',choices=list(ROOTS)+['all'])
    p.add_argument('--box',nargs=6,metavar=('X0','X1','Y0','Y1','B0','B1'))
    p.add_argument('--mesh',type=int,default=7)
    p.add_argument('--depth',type=int,default=25)
    p.add_argument('--limit',type=int,default=50000)
    p.add_argument('--single-alias',action='store_true',help='same alias evaluator as the Lean moderate certificate')
    p.add_argument('--export',type=Path,help='write exact parameter witnesses after successful replay')
    a=p.parse_args()
    if bool(a.root)==bool(a.box) or not 0<=a.mesh<=16:
        p.error('supply exactly one of --root and --box and a valid mesh exponent (0–16)')
    roots=list(ROOTS) if a.root=='all' else [a.root] if a.root else ['custom']
    witness={'format':'shallow-cross-helmholtz-v1','inputScaleBits':Q.bit_length()-1,
        'singleAlias':a.single_alias,
        'spectralMesh':{'numerators':[3*j for j in range((1<<a.mesh)+1)],'denominator':4*(1<<a.mesh)},
        'parameterAxes':['x','y','b=kappa^2'],
        'prefixEncoding':{'0':['lower-margin numerator'],'1':['axis','cut numerator']},
        'warning':'Leaf labels were re-evaluated by Python, not Lean theorems.', 'roots':[]}
    for name in roots:
        b=exact_box(a.box if a.box else ROOTS[name])
        trace=[] if a.export else None
        ok,n,k,d,lo,why,cell=cover(b,a.mesh,a.depth,a.limit,trace,a.single_alias)
        print(name,'CERTIFIED' if ok else 'UNCERTIFIED',f'visits={n} leaves={k} depth={d}',
            f'lower={lo}/{Q}' if lo is not None else 'no accepted margin',why,simple(cell),flush=True)
        if not ok:
            return 1
        if a.export:
            witness['roots'].append({'name':name,'ranges':[[i.lo,i.hi] for i in (b.x,b.y,b.b)],
                'visits':n,'leaves':k,'depth':d,'minLowerMarginNumerator':lo,'tree':trace})
    if a.export:
        a.export.write_text(json.dumps(witness,separators=(',',':'))+'\n')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
