#!/usr/bin/env python3
"""Additional exact continuous boxes for 16<=x<=24, 56<=y<=80.

Identical signed-mass formulas and outward-integer jets to
moderate_helmholtz.py. The only analytic change is retaining 64 rather than
24 terms in the nonnegative kappa^2 alias series. The scaled tail is
2 exp(128) (1/256)^64/(1-1/256) < 2^-100.
"""
import argparse
import json
from pathlib import Path
from moderate_helmholtz import I,J,Q,ZERO,ONE,Model as Base,exact_box,children,split_axis,simple,mesh


def calcH64(t,b):
    z=moment=cur=ONE
    accum=ZERO
    for j in range(1,128):
        z=z*t/j
        moment=moment+z
        if j%2:
            accum=accum+cur*moment
            cur=cur*b
    leading=2*accum+I(0,Q>>100)
    e=(-t).exp()
    den=(1-e).sq()
    return I(leading.lo,(leading*(1+2*e/den)).hi),3*e/den


class Model(Base):
    __slots__=()
    def __init__(self,x,y,b):
        if not (16*Q <= x.lo <= x.hi <= 24*Q and 56*Q <= y.lo <= y.hi <= 80*Q
                and 0<=b.lo<=b.hi<=Q//256):
            raise ValueError("outside 16..24 x 56..80")
        c=(y-x)/2
        m=(y+x)/2
        self.x,self.y,self.c,self.m=x,y,c,m
        self.xs,self.ms,self.ys=calcH64(x,b),calcH64(m,b),calcH64(y,b)
        self.em=(-c).exp()
        self.c2=c.sq()


def one(box,nodes):
    m=Model(box.x,box.y,box.b)
    try:
        return m.high()+sum(m.low(p)*p.span()//Q for p in nodes)
    except ArithmeticError:
        return None


def cover(box,nodes,cap,limit):
    todo=[(box,0)]
    n=k=h=0
    minval=None
    trace=[]
    while todo:
        b,d=todo.pop()
        n+=1
        h=max(h,d)
        val=one(b,nodes)
        if val is not None and val>0:
            k+=1
            trace.append([0,val])
            minval=val if minval is None else min(minval,val)
        else:
            if d>=cap or n>=limit:
                return False,n,k,h,minval,trace,simple(b)
            axis=split_axis(b)
            l,r=children(b,axis)
            trace.append([1,axis,(l.x,l.y,l.b)[axis].hi])
            todo.append((r,d+1))
            todo.append((l,d+1))
    return True,n,k,h,minval,trace,None


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--mesh',type=int,default=6)
    a.add_argument('--limit',type=int,default=40000)
    a.add_argument('--export',type=Path)
    args=a.parse_args()
    nodes=mesh(args.mesh)
    all_roots=[]
    for l,h in ((56,64),(64,80)):
        b=exact_box((16,24,l,h,0,'1/256'))
        ok,n,k,d,v,tree,why=cover(b,nodes,20,args.limit)
        print(l,h,'CERTIFIED' if ok else 'UNCERTIFIED',n,k,d,f"{v}/{Q}",why,flush=True)
        if not ok: return 1
        all_roots.append(dict(ranges=[[p.lo,p.hi] for p in (b.x,b.y,b.b)],tree=tree))
    if args.export:
        args.export.write_text(json.dumps(dict(inputScaleBits=112,cutoff='29/32',meshPow=args.mesh,roots=all_roots),separators=(',',':'))+'\n')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
