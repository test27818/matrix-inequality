#!/usr/bin/env python3
"""Original J-kernel moderate crossed rectangles: integer Darboux lower bounds.
Extends shallow_cross_helmholtz.py to 3<=x<=16, x<=y<=32, x<y.
"""
import argparse
from fractions import Fraction
import shallow_cross_helmholtz as s
from moderate_helmholtz import Q,I,exact_box, simple

class CompactModel(s.Model):
    def __init__(self,x,y,b,single_alias=False):
        c=(y-x)/2;m=(y+x)/2
        if not(3*Q<=x.lo<=x.hi<=16*Q and x.hi<y.lo<=y.hi<=32*Q and 0<=b.lo<=b.hi<=Q//256):
            raise ValueError('only separated cross boxes at depths 3..32')
        self.x,self.y,self.c,self.m=x,y,c,m
        hcol=(lambda t: (*s.calc_h(t,b),1)) if single_alias else (lambda t: s.hyperbolic_column(t,b))
        self.xs,self.ms,self.ys=map(hcol,(x,m,y))
        self.em=(-c).exp();self.c2=c.sq()
s.Model=CompactModel

def check(x0,x1,y0,y1,beta,mesh,depth,limit,b0='0',b1='1/256'):
    beta=Fraction(beta)
    if not Fraction(1,2)<=beta<=Fraction(29,32):raise ValueError('cutoff in [1/2,29/32]')
    s.BETA=I.q(beta)
    if s.BETA.lo!=s.BETA.hi or s.BETA.lo%(1<<mesh):raise ValueError('cutoff/mesh is not exact dyadic')
    box=exact_box((x0,x1,y0,y1,b0,b1))
    ok,n,k,d,lo,why,cell=s.cover(box,mesh,depth,limit)
    print('box',x0,x1,y0,y1,'b',b0,b1,'beta',str(beta),'mesh',mesh,
          'CERTIFIED' if ok else 'UNCERTIFIED',f'visits={n} leaves={k} depth={d}',
          f'lower={lo}/{Q}' if lo is not None else 'no accepted margin',why,simple(cell),flush=True)
    return ok
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('x0');p.add_argument('x1');p.add_argument('y0');p.add_argument('y1')
    p.add_argument('beta',nargs='?',default='3/4');p.add_argument('mesh',nargs='?',type=int,default=6)
    p.add_argument('--depth',type=int,default=25);p.add_argument('--limit',type=int,default=16000)
    p.add_argument('--b0',default='0');p.add_argument('--b1',default='1/256')
    a=p.parse_args()
    raise SystemExit(0 if check(a.x0,a.x1,a.y0,a.y1,a.beta,a.mesh,a.depth,a.limit,a.b0,a.b1) else 1)
