#!/usr/bin/env python3
"""Exact rational envelope comparisons for the finite-order outer branch."""
from fractions import Fraction as F
from math import factorial
q=F(27,10)
assert sum(F(1, factorial(j)) for j in range(6)) > q
assert sum(F(1,factorial(j)) for j in range(6))+F(1,720)/(1-F(1,7)) < F(14,5)
assert 54/q**7 < F(3,50)
assert F(19,2)*F(3,50) < 1
assert 178*F(3,50)+F(3,2)*18**2/q**8 < 20
assert F(14,5)**3 < F(47,10)**2
assert 3*F(9,8)*F(47,10) < 16
assert 3*F(26)**2/24 < 85
assert F(1,1)/(2*F(9,10)) < 1
assert F(20,1)/(2*F(9,10))+F(1,1)/(4*F(9,10)**3) < 12
assert 12/F(19,10)**2+2/F(19,10)**3 < 4
assert 85+2*16+4*16 < 192
assert 3*18/q**15 < 1
assert sum(F(1,factorial(2*j+1)) for j in range(4)) + F(1,factorial(9))/(1-F(1,110)) < F(59,50)
assert sum(F(1,factorial(2*j)) for j in range(4)) + F(1,factorial(8))/(1-F(1,90)) < F(31,20)
assert 1/(1-q**-16)**2 < F(101,100)
assert F(101,100)*(1-F(50,59)+F(11,20)/16+(1+F(31,20)/16)/255) < F(1,4)
assert 1+F(24,28) == F(13,7)
assert (2+4*F(13,7)+F(24,28))/(16*q**2) < F(1,10)
assert 2*(F(1,4)+96)/(28*q**14) < F(1,1000)
# At the left endpoint y=80, beta=28.  exp(-99/4) is bounded by exp(-24).
m=F(80+28,2)
term_m=(3*m*(1+F(80**2,8))+F(3,2)*m*m)/q**24
term_y=3*(80+2)*(F(3,8)*80**2+4*80+15)/q**31
assert (term_m+term_y)/(28*16) < F(1,1000)
assert 1-F(1,112)-F(192,28**2)-F(1,1000)-F(1,10)-F(1,1000) > F(1,2)
print("All outer moderate branch rational comparisons verified.")
