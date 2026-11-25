# Shallow finite-kernel outer branch: 3 ≤ x ≤ 8, y ≥ 48

This branch concerns the actual finite Hermite coefficients with 0≤κ≤1/16, including physical values κ=1/N. It does not substitute the κ=0 limit. The notation T,G,E,W, c=(y−x)/2, and m=(x+y)/2 follows [the other outer branch](outer-moderate.md).

## Bounded part 48 ≤ y ≤ 80

The four rectangles obtained by splitting [3,8] at 5 and [48,80] at 64 all pass [the exact integer checker](../../certificates/far_shallow_bridge.py), uniform in κ²∈[0,1/256]. For the shallow x column retain three aliases as justified in [shallow cross reconstruction](shallow-cross-helmholtz.md); for the deeper y and midpoint columns use 64 odd terms of the first alias as justified in [the other outer branch](outer-moderate.md). The cutoff is β=3/4, so the high envelope applies (βm≥1). The signed-mass identity is unchanged. [Exact replay：从源码重放](../../replay.sh) has 40, 40, 8, and 8 accepted leaves, respectively, with strictly positive lower bounds. All boxes are continuous and closed.

## Unbounded part y ≥ 80

Here c≥36. As in the other outer proof, write D=e^{-cp}G+E on 0≤p≤1/2, where G=T_x/(1+√W_x) and E=-T_m+Ze^{cp}T_y, Z=√W_x/(1+√W_y). Put g₀=G(0)≥xe^{-x}, by the same integral proof as §2.1 of [the other outer branch](outer-moderate.md). For v≥16 all the W and square-root bounds from that proof remain valid. Below we prove fresh constants for x in [3,8]. All numerical comparisons can be rerun in [the exact rational script](../../certificates/verify_far_shallow_bridge_constants.py).

### Shallow low-band derivatives

For positive p with the regularized value at zero let
$$
T_x(p)=(1-p)H_x+A_x(p)-B_x(p)
$$
where
$$
H_x=2\int_x^\infty\frac{e^{-v}\sinh(κv)/κ}{(1-e^{-v})^2}\,dv,\quad
A_x(p)=\int_x^\infty\frac{e^{-(2-p)v}}{(1-e^{-v})^2}\,dv,\quad
B_x(p)=\int_x^\infty\frac{(2-e^{-v})e^{-(1+p)v}}{(1-e^{-v})^2}\,dv.
$$
This is just the exact pole subtraction in the physical spectrum. Both derivatives of A are nonnegative, B′ is nonpositive and B″ is nonnegative, by differentiating exponentially convergent integrals. For x≥3, (1−e^{-x})^{-2}<28/25; for x≤8, κx≤1/2 and the exact integral of e^{-v} sinh(κv)/κ implies H_x≤(12/5)(x+6/5)e^{-x}. Indeed S(x)/x≤21/20, cosh(κx)≤113/100 and κ²≤1/256. Use
$$
\begin{split}
e^{px}A_x &\le(3/4)e^{-x},\
e^{px}A_x'&\le(3/4)(x+2/3)e^{-x},\
e^{px}|B_x'|&\le(56/25)(x+1)e^{-x},\
e^{px}\max(A_x'',B_x'')&\le(56/25)(x^2+2x+2)e^{-x}.
\end{split}
$$
These all follow by integrating e^{-av}v^j for j=0,1,2 and replacing the denominator by its x endpoint. Since A′−B′≥0, |T_x′|≤max(H_x,A_x′−B_x′), and |T_x″|≤max(A_x″,B_x″). After dividing by g₀≥xe^{-x}, and preserving the additional e^{-x/2} in the A and A′ bounds before multiplying by e^{px}, one obtains
$$
|T_x|,|T_x'|\le4g_0,\qquad |T_x''|\le25g_0.
$$
For T_x use T_x≥0 and T_x≤(1−p)H_x+A_x. This sign also gives W_x≤1.

Let (L(x)=(12/5)(x+6/5)e^{-x/2}).  On 0≤p≤1/2, the term p(1−p)e^{px} is increasing, so p(1−p)e^{px}H_x≤L(x)/4. Thus
$$
1-W_x\le L(3)/4+(3/8)e^{-3}<3/5,
$$
and W_x>2/5. Also e^{px}|T_x'|≤L(x), since A′−B′ times e^{px} is ≤3(x+1)e^{-x}≤L(x). Differentiating W gives
$$
|W_x'|\le(3/2+x/2)L(x)+(3/4)(1+x/2)e^{-x}<8.
$$
Both bounds decrease in x from 3. For its second derivative an upper bound is
$$
 (x^2/2+3x+2)L(x)
 +(3/4)(x^2/2+2x)e^{-x}+(28/25)(x^2+2x+2)e^{-x}<60.
$$
For the first polynomial times e^{-x/2}, bound the monomials of degrees 3 and 2 at their maxima 6 and 4, those of degrees 0 and 1 at 3. For the other two polynomials the maximum is at x=3, term by term. These are the exact bounds tested in the companion script.

Thus for u=√W_x we may take u≥3/5, |u′|≤7 and |u″|≤125. For C=(1+u)^{-1}, |C′|≤3 and |C″|≤74. Consequently |G″|≤(25+2·4·3+4·74)g₀<352g₀. At zero use the initial slope identity and bound as in the other outer proof, this time with (1−e^{-x})^{-2}≤28/25, S(x)/x≤21/20, cosh(κx)≤113/100 and x≥3. Then −T′(0)/T(0)<1/4, and T(0)<1 follows directly from H_x+A_x(0)<1. Thus |G′(0)|≤g₀/4 and |G′(1/2)|≤(1/4+176)g₀.

### High band and low-band correction

The same actual U,V are convex decreasing on [1/2,1−κ], with d/dp log J_v=−v−1/p+W_v′/W_v. At p=1/2 the shallow slope magnitudes obey |W_x′|/W_x≤20, and the mid and far ones obey |W_v′|/W_v<2. Thus −U′≤(m+4)U, −V′≤(m+13)V, where U,V≤2e^{-m/2}. The total high cost divided by c g₀ is therefore bounded by
$$
[2c/m+(4m+34)/c] e^{(x-c)/2}/x<1/100
$$
using m/c≤11/9, x≥3, x−c≤−28.

In the low E term, |Z|≤1, |Z′|≤8, |Z″|≤143 by combining shallow u bounds with the far inverse square-root factor's bounds. As in the other outer proof set w=(y+12)/2≥m+2. Using the uniform v≥8 bounds on T for m and y (DeepTail.lean), the cost numerator |E′(0)|+∫₀^{1/2}|c²E−E″|, times e^x, is at most
$$
[3w(1+y^2/8)+\tfrac32w^2] e^{(136-15y)/32}
+3(y+2)(y^2/4+11y/2+99)e^{8-11y/16}.
$$
For the second term use |(Ze^{cp}T_y)'|≤3(y+2)(c+9)e^{-11y/16} and the second derivative bound with bracket c²+18c+y+180; these follow by ordinary product rules from the derivative bounds just given. Use c≤y/2 to obtain the polynomial in the display. At y=80 the exponents are at most −33 and −47; all polynomial terms decrease when multiplied by their respective exponentials on y≥80. Divide by c x≥36·3 using g₀≥xe^{-x}. The low correction cost divided by c g₀ is <1/1000.

### Transform conclusion

The affine-extension identity proved in §2.3 of the other outer branch now gives a normalized transform lower bound
$$
1-1/144-352/36^2-1/1000-1/100-1/1000>2/3.
$$
The five deductions are exactly: the first slope, normalized weighted curvature, affine tail, high physical spectrum with its cutoff slope, and low E correction. Hence the true cosine transform is positive at every real frequency. This closes all 3≤x≤8, y≥48 when combined with the bounded integer proof. It does not cover x<3.
