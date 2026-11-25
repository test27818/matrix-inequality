#!/usr/bin/env python3
"""The Helmholtz evaluator with correlated (x,c,b), y=x+2c, m=x+c.

Reuses the single-alias exact arithmetic. Pass only dyadic cuts below 15/16.
Covers a fixed radial domain; a reject is not a counterexample.
"""
from fractions import Fraction
import argparse
from pathlib import Path
import json
import moderate_helmholtz as core

Q=core.Q

class Model(core.Model):
    def __init__(self,x,c,b):
        if not (16*Q <= x.lo <= x.hi <= 24*Q and 0 < c.lo <= c.hi <= 16*Q
                and 0 <= b.lo <= b.hi <= Q//256):
            raise ValueError('affine radial box outside supported domain')
        self.x,self.c,self.m,self.y=x,c,x+c,x+2*c
        self.xs,self.ms,self.ys=core.calc_h(x,b),core.calc_h(self.m,b),core.calc_h(self.y,b)
        self.em=(-c).exp()
        self.c2=c.sq()


def score(box,mesh):
    M=Model(box.x,box.y,box.b)
    try:
        return M.high()+sum(M.low(p)*p.span()//Q for p in mesh)
    except ArithmeticError:
        return None


def cover(box,mesh,dmax,limit,trace=None):
    todo=[(box,0)]
    visits=leaves=depth=0
    lo=None
    while todo:
        if visits>=limit:
            return False,visits,leaves,depth,lo,'limit',todo[-1][0]
        node,d=todo.pop()
        visits+=1
        depth=max(depth,d)
        v=score(node,mesh)
        if v is not None and v>0:
            leaves+=1
            lo=v if lo is None else min(v,lo)
            if trace is not None:
                trace.append([0,v])
            continue
        if d>=dmax:
            return False,visits,leaves,depth,lo,'depth',node
        axis=core.split_axis(node)
        l,r=core.children(node,axis)
        if trace is not None:
            trace.append([1,axis,(l.x,l.y,l.b)[axis].hi])
        todo.append((r,d+1))
        todo.append((l,d+1))
    return True,visits,leaves,depth,lo,'done',None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',choices=['left','right','all'],default='all')
    p.add_argument('--gap-min',type=Fraction,default=Fraction(4))
    p.add_argument('--gap-max',type=Fraction,default=Fraction(8))
    p.add_argument('--beta',type=Fraction,default=Fraction(29,32))
    p.add_argument('--mesh',type=int,default=6)
    p.add_argument('--depth',type=int,default=20)
    p.add_argument('--limit',type=int,default=20000)
    p.add_argument('--export',type=Path)
    a=p.parse_args()
    if not (Fraction(1,16)<=a.beta<Fraction(15,16) and 0<=a.mesh<=16 and
            Q % (a.beta.denominator*(1<<a.mesh))==0):
        p.error('β must be dyadic, at least 1/16 and strictly below 15/16, with exact 2^-112 spectral mesh')
    if not 0<a.gap_min<=a.gap_max<=16:
        p.error('require 0<gap-min<=gap-max<=16 (gap variable is c=(y−x)/2)')
    core.BETA=core.I.q(a.beta)
    mesh=core.mesh(a.mesh)
    data={'format':'compact-gap-helmholtz-v1','scaleBits':Q.bit_length()-1,
          'axes':['x','c=(y-x)/2','b=kappa^2'],'beta':str(a.beta),'meshDepth':a.mesh,
          'prefixEncoding':{'0':['scaled lower margin numerator'],'1':['axis','cut numerator']},
          'roots':[]}
    roots=[('left',16,20),('right',20,24)]
    for name,x0,x1 in roots:
        if a.root not in ('all',name):
            continue
        b=core.exact_box((x0,x1,a.gap_min,a.gap_max,0,'1/256'))
        trace=[] if a.export else None
        ok,n,k,d,lo,why,cell=cover(b,mesh,a.depth,a.limit,trace)
        print(name,'CERTIFIED' if ok else 'UNCERTIFIED',f'visits={n} leaves={k} depth={d}',
            f'lower={lo}/{Q}' if lo is not None else 'no accepted margin',why,core.simple(cell),flush=True)
        if not ok:
            return 1
        if a.export:
            data['roots'].append({'name':name,'ranges':[[v.lo,v.hi] for v in (b.x,b.y,b.b)],
                'tree':trace,'visits':n,'leaves':k,'depth':d,'minLowerMarginNumerator':lo})
    if a.export:
        a.export.write_text(json.dumps(data,separators=(',',':'))+'\n')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
