"""Cauchy-tail checks for compact_mid_spectrum.cpp (same 29/32 cutoff as 16-24 case)."""
from fractions import Fraction as F
from math import factorial
h=F(3,16) # upper bound, actual 3/32+1/16=5/32
def neg(x):
    # e^-x <= reciprocal of its positive Taylor lower sum.
    return 1/sum(x**i/factorial(i) for i in range(120))
for xmin,k in [(16,3),(12,4),(8,7),(6,9),(4,15),(3,23)]:
    t=F(xmin-1)
    hialias=4*33**3/(1-h)*neg((k-h)*t)/(1-neg(t))
    hidegree=2*k*F(11,4)**(33+k)*h**62/(1-h**2)
    loalias=(2*(1+F(1,16))*(33+F(1,k+1))+4*F(11,4)**3)*neg(k*t)/(1-neg(t))
    a=512*(hialias+hidegree); b=512*loalias
    assert max(a,b)<F(1,2**33)
    print('minimum x=%d aliases=%d: 2^-33 mixed-jet allowances verified' %(xmin,k))
