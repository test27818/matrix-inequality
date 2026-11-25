#!/usr/bin/env python3
"""Exact bounded bridge for x in [3,8] and y in [48,80]."""
from itertools import product
import argparse
import json
from pathlib import Path
from moderate_helmholtz import Q,I,exact_box,children,split_axis,simple
from shallow_cross_helmholtz import Model as ShallowModel,BETA,hyperbolic_column
from outer_moderate_helmholtz import calcH64


class Model(ShallowModel):
    __slots__=()
    def __init__(self,x,y,b):
        if not (3*Q<=x.lo<=x.hi<=8*Q and 48*Q<=y.lo<=y.hi<=80*Q and 0<=b.lo<=b.hi<=Q//256):
            raise ValueError("outside bridge")
        c=(y-x)/2
        m=(x+y)/2
        self.x,self.y,self.c,self.m=x,y,c,m
        self.xs=hyperbolic_column(x,b)
        self.ms=(*calcH64(m,b),1)
        self.ys=(*calcH64(y,b),1)
        self.em=(-c).exp()
        self.c2=c.sq()


def one(box,nodes):
    m=Model(box.x,box.y,box.b)
    try:
        return m.high()+sum(m.low(p)*p.span()//Q for p in nodes)
    except ArithmeticError:
        return None


def cover(box,ns):
    pending=[(box,0)]
    trace=[]
    n=leaves=depth=0
    best=None
    while pending:
        b,d=pending.pop()
        n+=1
        depth=max(depth,d)
        v=one(b,ns)
        if v is not None and v>0:
            leaves+=1
            best=v if best is None else min(v,best)
            trace.append([0,v])
        elif n>=20000 or d>=21:
            return False,n,leaves,depth,best,trace,simple(b)
        else:
            i=split_axis(b)
            l,r=children(b,i)
            trace.append([1,i,(l.x,l.y,l.b)[i].hi])
            pending.append((r,d+1))
            pending.append((l,d+1))
    return True,n,leaves,depth,best,trace,None


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--mesh',type=int,default=6)
    a.add_argument('--export',type=Path)
    args=a.parse_args()
    w=BETA.hi//(1<<args.mesh)
    nodes=[I(j*w,(j+1)*w) for j in range(1<<args.mesh)]
    roots=[]
    for (x1,x2),(y1,y2) in product(((3,5),(5,8)),((48,64),(64,80))):
        b=exact_box((x1,x2,y1,y2,0,'1/256'))
        ok,n,l,h,v,t,why=cover(b,nodes)
        print(x1,x2,y1,y2,'CERTIFIED' if ok else 'UNCERTIFIED',n,l,h,f'{v}/{Q}',why,flush=True)
        if not ok: return 1
        roots.append(dict(ranges=[[q.lo,q.hi] for q in (b.x,b.y,b.b)],tree=t))
    if args.export:
        args.export.write_text(json.dumps(dict(inputScaleBits=112,cutoff='3/4',meshPow=args.mesh,roots=roots),separators=(',',':'))+'\n')
    return 0


if __name__=='__main__': raise SystemExit(main())
