#!/usr/bin/env python3
"""Verify the scalar overestimates in deep-continuum-closure.md by rational bounds."""
from fractions import Fraction as F
from math import factorial
from moderate_helmholtz import I,Q
def up(name,val,cap):
    if val.hi>I.q(cap).lo: raise ArithmeticError(name)
    print('PASS',name)
e24=I.q(-24).exp()
up('small W/t^3 aliases: each Taylor lane',16*e24,I.q(F(1,10**7)))
up('small square numerator aliases',7*33*e24,I.q(F(1,10**6)))
up('large W aliases: each Taylor lane',26*e24,I.q(F(1,10**7)))
up('small E,C Taylor lane tails',I.q(2*F(3**32*34*33*2,factorial(36))),I.q(F(1,10**22)))
assert 2*F(9**84,factorial(85)) < F(1,10**40)
print('PASS main positive series tail')
up('linear alias sup on [0,1], including first two derivatives',1324*9*I.q(-15).exp()/(1-2*e24),I.q(F(1,250)))
up('glued L1 mass beyond p=1',e24,I.q(F(1,10**6)))
up('glued second variation beyond p=1',2*e24,I.q(F(1,10**6)))
up('boundary layer T>=96: mass and second derivative integrals',
  20_000_000*I.q(-78).exp()*(108**4+4*108**3+12*108**2+24*108+24),I.q(F(1,10**8)))
