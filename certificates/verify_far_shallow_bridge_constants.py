#!/usr/bin/env python3
"""Exact arithmetic for [3,8] × [80,∞) high-order varying-radius proof."""
from fractions import Fraction as F
from math import factorial
q=F(27,10)
s=F(1,2)
assert sum(F(1,factorial(j)) for j in range(6))>q
assert sum(F(3,2)**j/factorial(j) for j in range(8))>F(40,9)
assert sum(s**(2*j)/factorial(2*j+1) for j in range(4))+s**8/factorial(9)/(1-F(1,440))<F(21,20)
assert sum(s**(2*j)/factorial(2*j) for j in range(4))+s**8/factorial(8)/(1-F(1,360))<F(113,100)
assert (1-q**-3)**-2<F(28,25)
assert 2*F(28,25)*F(256,255)*F(21,20)<F(12,5)
assert 2*F(28,25)*F(256,255)*F(113,100)<F(72,25)
L3=F(12,5)*F(21,5)*F(9,40)
assert L3/4+F(3,8)/q**3<F(3,5)
assert 3*L3+F(3,4)*F(5,2)/q**3<8
first=F(12,5)*(F(108)/q**3 + F(18,5)*16/q**2 + F(28,5)*F(3,4)+F(12,5)/4)
last=(F(3,4)*F(21,2)+F(28,25)*17)/q**3
assert first+last<60
assert F(12,5)*F(7,5)+F(3,4)/12<4
assert F(3,4)*(1+F(2,9))/4+F(56,25)*F(4,3)<4
assert F(56,25)*F(32,3)<25
assert F(12,5)*F(21,5)/q**3+F(3,4)/q**6<1
assert F(28,25)*(1-F(20,21)+F(13,300)+(1+F(113,300))/255)<F(1,4)
assert F(8,1)/F(6,5)<7
assert F(60,1)/F(6,5)+F(64,1)/(4*F(3,5)**3)<125
assert 7/F(8,5)**2<3
assert 125/F(8,5)**2+98/F(8,5)**3<74
assert 25+2*4*3+4*74<352
assert (2+4*F(11,9)+F(34,36))/(3*q**14)<F(1,100)
assert 2*(F(1,4)+176)/(36*q**18)<F(1,1000)
m=46
tm=(3*m*(1+F(80**2,8))+F(3,2)*m*m)/q**33
ty=3*82*(F(80**2,4)+F(11,2)*80+99)/q**47
assert (tm+ty)/(36*3)<F(1,1000)
assert 1-F(1,144)-F(352,36**2)-F(1,1000)-F(1,100)-F(1,1000)>F(2,3)
print('All far-shallow bridge rational comparisons verified.')
