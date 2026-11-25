#!/usr/bin/env python3
"""Exact vertical section audit for the radial proof coverage table.

This checks only geometrical coverage, not any analytic bound or replay.
Each row is a closed x-slab and an intersection of affine bounds on y.
The intended domain is x>=0, y>=x.  x=0 can be ignored in the matrix
argument, but included if covered; the line y=x is covered separately.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations

@dataclass(frozen=True, order=True)
class Line:
    a: int
    b: int
    def __call__(self, x):
        return self.a*x+self.b

@dataclass(frozen=True)
class Row:
    name: str
    x0: int
    x1: int | None
    low: tuple[Line, ...]
    high: tuple[Line, ...]

DIAG=Line(1,0)
def box(name,x0,x1,y0,y1):
    return Row(name,x0,x1,(Line(0,y0),),() if y1 is None else (Line(0,y1),))
def gap(name,x0,x1,d0,d1):
    return Row(name,x0,x1,(Line(1,d0),),() if d1 is None else (Line(1,d1),))

PROVED=(
    Row("same-radius",0,None,(DIAG,),(DIAG,)),
    box("finite shallow square",0,4,0,4),
    box("zero-mixed complete shallow bridge (≤24)",0,3,0,24),
    box("zero-mixed and far, 24–48 plus ≥48",0,3,24,None),
    box("compact separated roots 0,1",3,8,16,24),
    box("compact 24–48 plus both far-shallow bridges",3,16,24,None),
    box("moderate plus outer-moderate",16,24,32,None),
    box("compact continuum low3",3,8,3,16),
    box("compact continuum mid8",8,12,8,16),
    box("compact finite phase c1",8,12,16,24),
    gap("compact continuum mid12",12,16,0,12),
    gap("continuum near full box x16–24",16,24,0,8),
    gap("correlated finite two roots",16,24,8,16),
    gap("correlated finite left extension",16,20,7,8),
    gap("deep continuum near strip",24,None,0,18),
    gap("deep separation",24,None,18,None),
)

def pieces(rows,x):
    work=[]
    for r in rows:
        if x<r.x0 or (r.x1 is not None and x>r.x1):
            continue
        lo=max((DIAG,*r.low),key=lambda l:l(x))
        hi=min(r.high,key=lambda l:l(x)) if r.high else None
        if hi is None or lo(x)<=hi(x):
            work.append((lo,hi))
    work.sort(key=lambda entry:entry[0](x))
    # Since the diagonal row is always present, every reported hole is open.
    gaps=[]
    end=DIAG
    for lo,hi in work:
        if end is None:
            break
        if lo(x)>end(x):
            gaps.append((end,lo))
        if hi is None or hi(x)>end(x):
            end=hi
    if end is not None:
        gaps.append((end,None))
    return tuple(gaps)

def arrangement(rows):
    cuts={F(0)}
    lines={DIAG}
    for r in rows:
        cuts.add(F(r.x0))
        if r.x1 is not None:
            cuts.add(F(r.x1))
        lines.update(r.low)
        lines.update(r.high)
    for u,v in combinations(lines,2):
        if u.a!=v.a:
            x=F(v.b-u.b,u.a-v.a)
            if x>=0:
                cuts.add(x)
    cuts=sorted(cuts)
    for i,s in enumerate(cuts):
        yield ("point",s,s,pieces(rows,s))
        t=cuts[i+1] if i+1<len(cuts) else None
        m=(s+t)/2 if t is not None else s+1
        yield ("open",s,t,pieces(rows,m))

def fmt(l):
    if l is None:
        return "∞"
    if l.a==0:
        return str(l.b)
    assert l.a==1
    return "x" if l.b==0 else f"x{l.b:+d}"

def summarize(rows):
    estate=[]
    for kind,a,b,g in arrangement(rows):
        if kind=="point":
            g=tuple((Line(0,l(a)),None if u is None else Line(0,u(a))) for l,u in g)
        estate.append((kind,a,b,g))
    for i,(kind,a,b,g) in enumerate(estate):
        if not g:
            continue
        # Omit point sections inherited from the formula on an adjacent open slab
        # only when their exact holes really agree.
        if kind=="point":
            neighboring=[estate[j][-1] for j in (i-1,i+1) if 0<=j<len(estate)]
            if any(tuple((Line(0,l(a)),None if u is None else Line(0,u(a))) for l,u in z)==g for z in neighboring):
                continue
            name=f"x={a}"
        else:
            name=f"{a}<x<{b if b is not None else '∞'}"
            if i>1 and estate[i-2][0]=="open" and estate[i-2][3]==g:
                continue
            j=i+2
            while j<len(estate) and estate[j][0]=="open" and estate[j][3]==g:
                b=estate[j][2]
                j+=2
            name=f"{a}<x<{b if b is not None else '∞'} (internal endpoints audited separately)"
        print(name+": "+"; ".join(f"{fmt(l)}<y<{fmt(u)}" for l,u in g))
    return not any(g for _,_,_,g in estate)

if __name__=="__main__":
    ok=summarize(PROVED)
    print("geometry fully closed:",ok)
    raise SystemExit(0 if ok else 1)
