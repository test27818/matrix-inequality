#!/usr/bin/env python3
"""Exact scalar constants for far-shallow-physical-proof.md; no third-party packages.

Only rational arithmetic is used, including rigorous Taylor enclosures of e^x.
This checks constants; statements about variable parameters are proved in the note.
"""
from fractions import Fraction as Q


def erange(x, n=80):
    x = Q(x)
    if x < 0:
        lo, hi = erange(-x, n)
        return 1/hi, 1/lo
    assert x < n+2
    s = t = Q(1)
    for j in range(1, n+1):
        t *= x / j
        s += t
    return s, s + t*x/(n+1)/(1-x/(n+2))


def upper(x): return erange(x)[1]
def lower(x): return erange(x)[0]


def check(name, *tests):
    for test in tests:
        if not test: raise ArithmeticError(name)
    print('PASS', name)


def main():
    q,h,k=Q('.00036'),Q('.019'),Q(1,256)
    check('Taylor endpoints', upper(-8)<q, upper(-8)<Q('.00034'), upper(-4)<h,
          lower(1)>Q(27,10), lower(Q(7,2))>23, Q(27,10)**14 > 10**6,
          (upper(1)-lower(-1))/2<Q('1.18'), (upper(1)+upper(-1))/2<Q('1.55'))
    A,B,E=Q('2.77'),Q('.02'),Q('2.01')
    check('A,B,E columns',2*(Q('1.18')+Q('1.55')/8)/((1-k)*(1-q)) < A,
          h/(1-q)<B, 2/(1-q)<E)
    T,V=Q('2.81'),Q(37)
    check('shallow T jet', A+B/8<T, E/8<T, A<T, (B+E)*Q(9,8)<T,
          E*Q(73,4)<V, T*16*q<1)
    # Every polynomial subexpression in v multiplying exp(-v/2) below has
    # nonnegative coefficients and its monomials of degrees 1,2,3 decrease at v>=8.
    v=Q(8)
    z0,z1,z2=Q('.22'),Q('2.4'),Q(28)
    check('shallow Z jet', T/2*v*h < z0, T*(1+v/2+Q(1,2))*v*h < z1,
          (T*(2*v+v*v/2+2+v)+V/2)*v*h < z2)
    r=Q('.85'); L0,L1,L2=Q('.55'),Q('.18'),Q('.25')
    check('inverse-root jet',r*r<1-z0,1/(1+r)<L0,1/(2*r*(1+r)**2)<L1,
          1/(4*r**3*(1+r)**2)+1/(2*r**2*(1+r)**3)<L2)
    check('G curvature',V*L0+2*T*L1*z1+T*(L1*z2+L2*z1*z1)<42)
    check('G endpoint slope',
          2*(Q(9,59)+Q(11,160)+Q(9,8*256))/(1-k)<(1-q)/2)
    check('all-radius T jet',2/((1-k)*(1-q))<Q('2.04'),
          Q('2.04')+B/9<3, E/9<3, (B+E)<3,
          E<3, 1/(1-q)<Q('1.01'), 2/(1-q)<E)
    check('half comparison', Q('2.04')/2+Q('1.01')/9<Q('1.2'),
          Q('1.01')+E<Q('3.1'),
          Q('3.1')*9*h+Q('1.2')*126/23<8)
    v=Q(32); d=Q(1,10**6)
    check('deep Z jet',3*v*d<z0,6*v*(Q('1.5')+v/2)*d<1,
          (6*v*(2*v+v*v/2)+6*v*(2+v)+6*v*v)*d<1)
    check('polarization jets', z1/(2*r)<2, z2/(2*r)+z1**2/(4*r**3)<22,
          L1+L2<1,22+2*2+1<=28)
    # For Y>=50, after dividing by indicated powers of Y.
    Y=Q(50)
    check('Q jets', Q(3,2)+3/Y<=2,Q(3,4)+3/Y+3/Y<=2)
    c1,c2,c3=700,400,22000
    check('shallow exponential ratios',
          upper(Q(17,2))/8 < c1,upper(8)/8<c2,upper(12)/8<c3)
    a=Q(24); A=a+10
    e0,e1,e2=[Q(10,27)**i for i in (22,24,33)]
    r0=3*A*c1*e0+3*Y*c3*e2
    r1=3*A*c1*e0+(9*Y+2*Y**2)*c3*e2
    r2=3*A**2*c2*e1+(84*Y+12*Y**2+2*Y**3)*c3*e2
    err=a*r0/2+3*r1/(2*a)+3*r2/(4*a)
    check('unbounded remainder',Q(3,24)<Q(15,16),err<Q(3,10000))
    affine=2*(Q(1,4)+Q(45,2))*Q('.00034')/16
    H=(Q(9)+Q(24,16))/16
    margin=1-Q(1,64)-Q(45,256)-Q(1,1000)-Q(3,10000)-H
    check('affine, high and final budget',affine<Q(1,1000),H==Q(21,32),margin>Q(1,8))
    print('remainder cap =',err)
    print('normalized guaranteed margin =',margin,'> 1/8')

if __name__=='__main__': main()
