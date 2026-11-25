#!/usr/bin/env python3
"""Rational audit for deep-separation-proof.md (all x>=24, k*x<=1, c>=9)."""
from fractions import Fraction as Q
from verify_far_shallow_physical import upper, lower, check


def tail(n,b=Q(5,6),t=6):
    """Upper bound for int_t^inf u^n exp(-bu) du (IBP n times)."""
    if n==0: return upper(-b*t)/b
    return t**n*upper(-b*t)/b+n*tail(n-1,b,t)/b


def main():
    t=Q(6); v=Q(24); beta=Q(3,4); c=Q(9)
    zi=upper(-24); ex6=upper(-6); ex5=upper(-5)
    AA,BB,EE=map(Q,('2.04','1.01','2.01'))
    check('positive column caps', 2/((1-Q(1,24**2))*(1-zi))<AA,
          1/(1-zi)<BB,2/(1-zi)<EE,2*(Q('1.18')+Q('1.55')/24)/((1-Q(1,24**2))*(1-zi))<Q('2.77'))
    # T/(e^{k v}q) <= 3(t+1); |T'| corresponding <= 3.2v, |T''| <=4v².
    check('universal T jets', AA*Q(25,24)<3, BB<3, AA*Q(25,24)<Q('3.2'),
          (BB+EE)*Q(25,24)<Q('3.2'),(BB+EE)*(1+2/v+2/v**2)<4,
          1+2/v+2/v**2<Q('1.1'))
    check('W and Z kernel envelopes',3*(t+1)*ex5<Q('.15'),Q('.9')**2<Q('.85'),
          3*Q(25,24)<=4,3*Q(25,24)+Q('3.2')<=8,
          3*(1+2/v)<=5,3*(1+2/v)+(2+2/v)*Q('3.2')+4<=20)
    b=Q(5,6); b2=2*b
    B=4*(tail(1)+2*tail(0))
    C=5*(tail(2)+4*tail(1))
    B2=16*(tail(3,b2)+4*tail(2,b2)+4*tail(1,b2))
    check('integrated envelope B, mu C and mu B^2',B<Q(1,3),C<4,B2<Q(1,2))
    # First-order column sign: max(A, B'+abs(E')), not an unsigned sum.
    check('shallow true derivatives at x>=24',Q('2.77')+BB*ex6/6<Q('2.78'),
          (EE+BB*ex6)*Q(25,24)<Q('2.81'),
          Q('1.1')*(EE+BB*ex6)<Q('2.3'))
    check('uniform G second variation',
          Q('.55')*Q('2.3')+2*Q('.18')*Q('2.81')/3+
          Q('2.78')*(Q('.18')*4+Q('.25')/2)<4)
    check('zero derivatives and cutoff ordinate slope',
          3*(v+1)*upper(-23)<1,
          (1+1/(v*beta)+Q('3.2')*ex5)/Q('.85')<Q('1.3'))
    check('D error rates',Q('1.1')*(BB+EE)<Q('3.4'),
          2*Q('3.4')+Q('3.5')*10+Q(4,3)*Q('3.2')<47,
          2*Q('3.2')+Q(4,3)*Q('3.5')<12,
          47+Q('9.9')<57,12+Q('3.5')<16,7+Q('3.5')<11)
    # For c>=9, 1+c/12 <= (c/9)(1+9/12), etc.
    ecap=(1+c/12)*(57/c+16+11*c)*upper(1-Q(11,12)*c)
    high=(1/(beta**2*v)+Q('2.6')/beta*(1/c+1/v))*upper(6-beta*c)
    aff=Q('8.5')*upper(-beta*c)/c
    check('uniform half-gap costs',Q(2,9)<Q(11,12),ecap<Q('.151'),
          high<Q('.29'),aff<Q('.002'))
    margin=1-1/(4*c)-4/c-Q('.151')-Q('.29')-Q('.002')
    check('positive final physical margin',margin>Q(1,13))
    print('rational remaining margin =',margin,'> 1/13')
if __name__=='__main__': main()
