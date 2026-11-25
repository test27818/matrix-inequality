#!/usr/bin/env python3
"""Exact outward checks for the singular-endpoint window [0,3] × [24,48].

Uses an analytic subtraction on t<=3 plus independent radial moments.
Floats do not enter acceptance or subdivision. All tails are analytic.
"""
import argparse
from fractions import Fraction
from math import factorial
from moderate_helmholtz import Q, I, J, ZERO, ONE, exact_box, children, simple, calc_h, tjet, ceildiv
from shallow_cross_helmholtz import hyperbolic_column, expanded_tjet

A = I.q(3)
BETA = I.q(Fraction(1,2))
NF = 14
NE = 23
def radius_check():
    """Exact analytic tails. The first omitted term can pay its own tail
    with the factor 2 because every next-term ratio is at most 1/2."""
    j = NF+1
    if not (Fraction(9,(2*j+1)*(2*j+2)) < Fraction(1,2) and
            Fraction(36,(2*j+2)*(2*j+3)) < Fraction(1,2) and
            Fraction(3,2*(NE+2)) < Fraction(1,2)):
        raise ArithmeticError('unsupported truncation tails')
    den_tail = I.q(Fraction(2*3**(2*j-2),factorial(2*j)))
    f_tail = I.q(Fraction(12*(2*j+2)**2*3**(2*j-1),factorial(2*j+1)))
    e_tail = I.q(9*Fraction(3,2)**(NE+1)/factorial(NE+1))
    return I(0,den_tail.hi),I(-f_tail.hi,f_tail.hi),I(-e_tail.hi,e_tail.hi)

DEN_TAIL,F_TAIL,E_TAIL=radius_check()


class SmallRadial:
    __slots__=('t','b','h3','e','e3','mom_f','mom_e')
    def __init__(self,t,b,radial_power):
        self.t,self.b=t,b
        self.h3=hyperbolic_column(A,b)
        self.e=t.exp()
        self.e3=(t-A).exp()
        K=1<<radial_power
        width=(A-t)/K
        mom_f=[ZERO for _ in range(NF)]
        for j in range(K):
            v=I(((K-j)*t.lo+j*3*Q)//K,ceildiv((K-j-1)*t.hi+(j+1)*3*Q,K))
            v2=v.sq()
            power=ONE
            den=ZERO
            nums=[]
            for n in range(1,NF+1):
                den=den+power/factorial(2*n)
                nums.append((v*power)/factorial(2*n+1))
                power=power*v2
            den=den+DEN_TAIL
            for i,num in enumerate(nums):
                mom_f[i]=mom_f[i]+num/den
        self.mom_f=[x*width for x in mom_f]
        self.mom_e=[]
        tpow=ONE
        for n in range(1,NE+4):
            tpow=tpow*t
            self.mom_e.append((I.q(3**n)-tpow)/n)

    def eval(self,p):
        mu=1-p
        mu2=mu.sq()
        mu_old=b_old=ONE
        f0=f1=f2=ZERO
        for j,rad in enumerate(self.mom_f,1):
            mu_new=mu_old*mu2
            b_new=b_old*self.b
            f0=f0+(mu*(mu_new-b_new))*rad
            f1=f1+(-(2*j+1)*mu_new+b_new)*rad
            f2=f2+((2*j+1)*j*mu*mu_old)*rad
            mu_old,b_old=mu_new,b_new
        F=J(f0+F_TAIL,f1+F_TAIL,f2+F_TAIL)
        ek=ONE
        e0=e1=e2=ZERO
        for k in range(NE+1):
            e0=e0+ek*self.mom_e[k]
            e1=e1-ek*self.mom_e[k+1]
            e2=e2+(ek*self.mom_e[k+2])/2
            ek=(ek*(-p))/(k+1)
        E=J(e0+E_TAIL,e1+E_TAIL,e2+E_TAIL)
        return J(self.e3)*expanded_tjet(A,p,self.h3)+J(self.e)*(E-F)


class Model:
    __slots__=('x','y','c','m','b','rad','ms','ys','em','c2')
    def __init__(self,x,y,b,radial_power):
        if not (0 <= x.lo <= x.hi <= 3*Q and 24*Q <= y.lo <= y.hi <= 48*Q
                and 0 <= b.lo <= b.hi <= Q//256):
            raise ValueError('outside [0,3]×[24,48]×[0,1/256]')
        c=(y-x)/2
        m=(x+y)/2
        self.x,self.y,self.c,self.m,self.b=x,y,c,m,b
        self.rad=SmallRadial(x,b,radial_power)
        self.ms,self.ys=calc_h(m,b),calc_h(y,b)
        self.em=(-c).exp()
        self.c2=c.sq()

    @staticmethod
    def W(t,p,T):
        return 1-p*((p-1)*J(t)).e()*T

    def low(self, pp):
        x,m,y,c=self.x,self.m,self.y,self.c
        p=J(pp,ONE)
        tx,tm,ty=self.rad.eval(pp),tjet(m,pp,self.ms),tjet(y,pp,self.ys)
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
        Wx=self.W(x,p,self.rad.eval(BETA))
        Wy=self.W(y,p,tjet(y,BETA,self.ys))
        ep=(J(x*(1-BETA))-p*J(c)-J(x)*(p-J(BETA))).e()
        V=(ep/p)*(Wx*Wy).s()
        return (2*V.d-self.c2*ep.f/(m*BETA)).lo


def one(box,nodes,radial_power):
    try:
        M=Model(box.x,box.y,box.b,radial_power)
        return M.high()+sum(M.low(p)*p.span()//Q for p in nodes)
    except ArithmeticError:
        return None


def cover(box,power,maxdepth,limit,radial_power):
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
        v=one(cell,nodes,radial_power)
        if v is not None and v>0:
            accepted+=1
            lo=v if lo is None else min(lo,v)
            continue
        if d>=maxdepth:
            return False,visits,accepted,depth,lo,'depth',cell
        axis=max(range(3),key=lambda i:(16*cell.x.span(),2*cell.y.span(),320*cell.b.span())[i])
        L,R=children(cell,axis)
        todo.append((R,d+1))
        todo.append((L,d+1))
    return True,visits,accepted,depth,lo,'done',None


ROOTS={f'x{i}-y{j}':(i,i+1,lo,hi,0,'1/256') for i in range(3) for j,(lo,hi) in enumerate(((24,32),(32,48)))}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',choices=list(ROOTS)+['all'])
    p.add_argument('--box',nargs=6,metavar=('X0','X1','Y0','Y1','B0','B1'))
    p.add_argument('--mesh',type=int,default=6)
    p.add_argument('--radial-mesh',type=int,default=5)
    p.add_argument('--depth',type=int,default=25)
    p.add_argument('--limit',type=int,default=50000)
    a=p.parse_args()
    if (bool(a.root)==bool(a.box) or not 0<=a.mesh<=16 or not 0<=a.radial_mesh<=16):
        p.error('supply one of --root and --box and positive integer mesh exponents at most 16')
    roots=list(ROOTS) if a.root=='all' else [a.root] if a.root else ['custom']
    for name in roots:
        b=exact_box(a.box if a.box else ROOTS[name])
        ok,n,k,d,lo,why,cell=cover(b,a.mesh,a.depth,a.limit,a.radial_mesh)
        print(name,'CERTIFIED' if ok else 'UNCERTIFIED',f'visits={n} leaves={k} depth={d}',
            f'lower={lo}/{Q}' if lo is not None else 'no accepted margin',why,simple(cell),flush=True)
        if not ok:
            return 1
    return 0

if __name__=='__main__':
    raise SystemExit(main())
