#!/usr/bin/env python3
"""Exact signed-mass experiments confined to separated compact rectangles.

No certificate is asserted on boxes for which the checker stops unresolved.
Uses the physical moderate/shallow one-shallower-end proof, β=3/4, 112-bit
outward integer enclosures, and up to three geometric aliases at small depth.
"""
import argparse
from fractions import Fraction
import shallow_cross_helmholtz as shallow
import json
from pathlib import Path
from moderate_helmholtz import Q,I,exact_box,children,split_axis,simple
from shallow_cross_helmholtz import Model as ShallowModel,BETA,hyperbolic_column


class Model(ShallowModel):
    __slots__=()
    def __init__(self,x,y,b):
        assert 3*Q <= x.lo <= x.hi <= 24*Q
        assert x.hi < y.lo <= y.hi <= 32*Q
        assert 0 <= b.lo <= b.hi <= Q//256
        c=(y-x)/2
        m=(x+y)/2
        self.x,self.y,self.c,self.m=x,y,c,m
        self.xs=hyperbolic_column(x,b)
        self.ms=hyperbolic_column(m,b)
        self.ys=hyperbolic_column(y,b)
        self.em=(-c).exp()
        self.c2=c.sq()


def one(b,ns):
    m=Model(b.x,b.y,b.b)
    try: return m.high() + sum(m.low(p)*p.span()//Q for p in ns)
    except ArithmeticError: return None


def cover(root,ns,limit,cap):
    todo=[(root,0)]
    n=k=depth=0
    lower=None
    tr=[]
    while todo:
        b,d=todo.pop()
        n+=1
        depth=max(depth,d)
        v=one(b,ns)
        if v is not None and v>0:
            tr.append([0,v])
            k+=1
            lower=v if lower is None else min(lower,v)
        elif n>=limit or d>=cap:
            return False,n,k,depth,lower,tr,simple(b)
        else:
            i=split_axis(b)
            l,r=children(b,i)
            tr.append([1,i,(l.x,l.y,l.b)[i].hi])
            todo.append((r,d+1))
            todo.append((l,d+1))
    return True,n,k,depth,lower,tr,None


ROOTS=((3,5,16,24),(5,8,16,24),(8,12,20,24),(12,16,22,24),
       (16,18,26,32),(18,20,28,32),(20,22,30,32),(22,24,31,32),
       (3,4,9,16),(4,5,10,16),(5,6,11,16),(6,7,12,16),(7,8,13,16),
       (3,5,11,16),(5,8,14,16),(8,12,22,24))

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--mesh',type=int,default=6)
    a.add_argument('--limit',type=int,default=8000)
    a.add_argument('--depth',type=int,default=16)
    a.add_argument('--root',type=int)
    a.add_argument('--beta',default='3/4')
    a.add_argument('--export',type=Path)
    args=a.parse_args()
    beta=Fraction(args.beta)
    assert Fraction(1,2)<=beta<=Fraction(29,32)
    shallow.BETA=I.q(beta)
    assert shallow.BETA.hi==shallow.BETA.lo and shallow.BETA.hi % (1<<args.mesh)==0
    w=shallow.BETA.hi//(1<<args.mesh)
    ns=[I(j*w,(j+1)*w) for j in range(1<<args.mesh)]
    roots=[]
    for i,rr in enumerate(ROOTS):
        if args.root is not None and i!=args.root: continue
        b=exact_box((*rr,0,'1/256'))
        ok,n,k,h,v,tr,why=cover(b,ns,args.limit,args.depth)
        print(i,*rr,'CERTIFIED' if ok else 'UNCERTIFIED',n,k,h,f'{v}/{Q}',why,flush=True)
        if ok:
            roots.append(dict(box=rr,ranges=[[p.lo,p.hi] for p in (b.x,b.y,b.b)],tree=tr))
            if args.export: args.export.write_text(json.dumps(dict(inputScaleBits=112,cutoff=str(beta),meshPow=args.mesh,roots=roots),separators=(',',':'))+'\n')
    return 0
if __name__=='__main__': raise SystemExit(main())
