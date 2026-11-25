"""Exact rational majorants used in the compact continuum weighted-Peano chart."""
from math import factorial
from fractions import Fraction as F
# e^t lower and upper by positive Taylor, where the tail is bounded by a geometric series.
def el(z,n=110): return sum(F(z)**k/factorial(k) for k in range(n+1))
def eu(z,n=110): return el(z,n)+F(z)**(n+1)/factorial(n+1)/(1-F(z,n+2))
assert eu(15)<4_000_000 and el(3)>20 and eu(3)<21
# fixed real v<=24, |q'-q|<=1/8. On high-p cells |q'|<=5/8.
# P_d=(E_{d-2}(jv)+(jv)^{d-1} C_{d-1}+(jv)^d C_d)/(1-(q'/j)^2).
assert 2*(25+4_000_000*(24**2+24**3))<10**12
# Cauchy q Taylor coefficient order <=2 multiplies by <=64.
assert 64*8*10**12<10**16
# geometric indices t>=0; (n+3+t)^2 <= (n+3)^2 (1+t)^2.
# sum_{k>=0}(1+k)^2/20^k<2.
z=F(1,20); assert (1+z)/(1-z)**3<2
assert F(2*10**16*38**2,1)/el(105,280)<F(1,2**45)
# C_k(qv) tail at l=35+, qv<=12,v<=24, after up to two q differentiations.
assert 2*24**2*F(12**68,factorial(68))<F(1,2**45)
print("PASS: compact Peano alias, centered angular Taylor and derivative tails < 2^-45")
