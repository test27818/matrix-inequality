"""Exact coarse rational estimates for analytic tails in compact_near_spectrum.cpp."""
from fractions import Fraction as F
from math import factorial
# e < 11/4, e^15 > 3,000,000 (positive Taylor lower sum).
lo=sum(F(15)**i/factorial(i) for i in range(80))
assert lo > 3_000_000
# q can be at most 1/8 on high branch and has radius 1/16 in Cauchy.
qh=F(3,16)
# 4 * 33^3 * e^(3/16) * e^(-(3-3/16)*15)/(1-e^-15)
# e^(3/16) <= 1/(1-3/16), e^(45*3/16?) use e^3 as 45/16 < 3.
hi_alias=4*33**3/(1-qh)*F(11,4)**3*F(1,3_000_000)**3/(1-F(1,3_000_000))
# e^((1-j)*Re(t))*sum |jt|^i/i! <= e^36 for j<=3.
hi_degree=6*F(11,4)**36*qh**62/(1-qh**2)
# Low T aliases: all exponents differ from -(j-1)t by a factor <= e^(33/16)<(11/4)^3; 4 such bounds (2 plus 2).
lo_alias=(2*(1+F(1,16))*(33+F(1,4))+4*F(11,4)**3)*F(1,3_000_000)**3/(1-F(1,3_000_000))
# Each (t^i p^j)-Taylor coefficient for i,j<=2 times i! costs <=2*16^2.
print('scaled high- and low-tail dyadic ratios (must be <1):')
print(512*2**33*(hi_alias+hi_degree),512*2**33*lo_alias)
assert 512*(hi_alias+hi_degree)<F(1,2**33) and 512*lo_alias<F(1,2**33)
