# Original finite Hermite kernel: outer branch 16 ≤ x ≤ 24, y ≥ 56

This recovers the moderate outer strip for the **finite** kernel, not only its κ=0 limit. For the unbounded part assume x≤N, which follows in the lift, and 0≤κ=1/N≤1/16. In the bounded part the integer checker actually covers the whole continuous κ² interval [0,1/256].

## 1. Bounded part 56 ≤ y ≤ 80

Use the exact signed-mass proof in [moderate Helmholtz](moderate-helmholtz.md), with identical jets, high cut at 29/32 and all outward arithmetic unchanged. Extend its positive Taylor column from 24 to 64 terms: for 0≤t≤80 (even 128 is allowed here), the omitted q=1 alias is bounded by
$$
2e^{128}(1/256)^{64}/(1-1/256)<2^{-100}.
$$
Indeed e<3 and 3¹²⁸<2²⁰⁴. All moving-alias enclosures are uniform in t>0. The two complete closed parameter boxes [16,24]×[56,64] and [16,24]×[64,80] are passed by [the integer program](../../certificates/outer_moderate_helmholtz.py), with the [output：从源码重放](../../replay.sh) and split trees. Replay from the repository root:
```sh
python3 certificates/outer_moderate_helmholtz.py --mesh 6
```
Subdivision never substitutes center samples; the proof divides the whole box along closed coordinate hyperplanes and sums rigorous lower endpoints on each continuous spectral cell. The spectral hypothesis gives the actual mixed Hermite kernel by [Fourier sampling](../continuum/positive-periodization.md).

## 2. Unbounded part y ≥ 80

Use κ=1/N, c=(y−x)/2 ≥28, m=(x+y)/2=x+c, and the original actual profile D=J_m−√(J_x J_y). Write
$$
T_v(p)=e^{-pv}/p-J_v(p),\qquad W_v(p)=1-p e^{pv}T_v(p),
$$
where T_v has its regularized value at p=0. On 0≤p≤1/2 use the exact identity
$$
D=e^{-cp}G+E,\quad G=T_x/(1+\sqrt{W_x}),\quad
E=-T_m+\frac{\sqrt{W_x}}{1+\sqrt{W_y}}e^{cp}T_y.
$$
These identities and actual derivative conventions are also verified in ZeroMode.lean and LowJet.lean. Let g_0=G(0)>0.

**Real low- and high-frequency inputs.** For v≥8 the following paper bounds follow by differentiating the positive exponential alias series, with the same quantities checked by ShallowDerivative.lean, ZeroMode.lean and DeepTail.lean:
$$
|T_v|,|T_v'|\le3(v+2)e^{-15v/16},\qquad |T_v''|\le3(v+2)^2e^{-v}
\quad(0\le p\le1/2).
$$
For clarity, the 0th, first and second series expressions are T0, T1 and T2 respectively (not the half-second Taylor lane). The proof of the bound is in ShallowSeries.lean and DeepTail.lean, which stay in this archive and require no recovered C++.
For 1/2≤p≤1−κ let U=J_m and V=√(J_xJ_y). Both are nonnegative, convex and nonincreasing by the positive convolution method (PositiveConvolution.lean, HighSpectrum.lean, HighProfiles.lean). The original density bound J_v≤e^{-pv}/p is checked for every p>0 in SharpOrdinate.lean. The true first slopes satisfy
$$
-U'(1/2)-V'(1/2)\le (4m+24)e^{-m/2}
$$
by HalfSlope.lean and HalfProfiles.lean, valid for all x,y≥8. The cutoff derivative needed below is the true one-sided derivative; its explicit physical formula is independently checked in CutoffSlope.lean. Alternatively all integrations can first stop at r<1−κ, and the slope term has the stated variation bound independently of r.

### 2.1 Low-frequency main term

Put S(v)=sinh(κv)/κ, with S(v)=v for κ=0. The exact regularized values follow directly from the density after subtracting its e^{-pv} pole:
$$
T_x(0)=\int_x^\infty\frac{2e^{-v}(S(v)-1+e^{-v})}{(1-e^{-v})^2}\,dv,\qquad
T_x'(0)=-\int_x^\infty\frac{2e^{-v}(S(v)-v)}{(1-e^{-v})^2}\,dv.
$$
The first integrand is at least 2e^{-v}(S(x)+v−x−1) since S'≥1 and x≥16. Thus g₀=T_x(0)/2≥e^{-x}S(x)≥xe^{-x}. Also T_x(0)<1 by the bound above. At p=0, W_x=1 and W'_x=−T_x(0), giving
$$
G'(0)/g_0=T_x(0)/4+T_x'(0)/T_x(0).
$$
Integrate e^{-v}S(v) to obtain e^{-x}(S(x)+cosh(κx))/(1−κ²). Since κx≤1, S(x)/x≤sinh(1)<59/50 and cosh(κx)≤cosh(1)<31/20. Therefore
$$
0\le -T_x'(0)/T_x(0)\le
\frac{101}{100}\left(1-\frac{50}{59}+\frac{11}{20\cdot16}
+\frac{1+31/(20\cdot16)}{255}\right)<\frac14.
$$
The elementary estimates for sinh(1), cosh(1) follow by their positive Taylor series with next-term geometric tails. Consequently |G'(0)|≤g₀/4.

For any v≥16 and p≤1/2, use the three low-frequency bounds on T to get
$$
|e^{pv}T_v|\le 3(v+2)e^{-7v/16}<3/50,\quad
0.97<W_v\le1,\quad |W_v'|\le1,\quad |W_v''|\le20.
$$
Indeed the first factor decreases for v≥16. For the first derivative, multiply the same factor by 1+v/2+1/2, which still decreases after multiplying by the stated exponential. For the second derivative, use the sum
$$
3(v+2)(v^2/2+3v+2)e^{-7v/16}+\tfrac32(v+2)^2e^{-v/2}.
$$
Each monomial in either polynomial decreases against its exponential starting at v=16, so evaluation at 16 is sufficient. For W_v≤1 use T_v≥0, which follows from J_v≤e^{-pv}/p and continuity at zero. If u=√W_v and A=(1+u)^{-1}, then u≥9/10 and the ordinary square-root and reciprocal rules give |u'|≤1, |u''|≤12, |A'|≤1, |A''|≤4. For x≤24, the bounds with g₀≥xe^{-x} give
$$
\max(|T_x|,|T_x'|)/g_0\le16,\quad |T_x''|/g_0\le85,\quad
|G''|/g_0\le85+2\cdot16+4\cdot16<192.
$$
In particular |G'(1/2)|≤(1/4+96)g₀.

### 2.2 Unbounded upper radius and error budget

Let Z=√W_x/(1+√W_y). The preceding derivative rules show |Z|≤1, |Z'|≤2, |Z''|≤18. For s=m and w=(y+28)/2≥m+2, the contributions of −T_m to
( |E'(0)|+\int_0^{1/2}|c^2E-E''|\)
after multiplication by e^x are bounded by
$$
\big[3w(1+y^2/8)+\tfrac32 w^2\big] e^{(408-15y)/32}.
$$
This uses c≤y/2 and x−15m/16=(17x−15y)/32; replacing e^{-m} by e^{-15m/16} only enlarges the second derivative cost. For the other contribution, y≥80 and c≤y/2 imply e^{cp}e^{-15y/16}≤e^{-11y/16}. Product differentiation with the above bounds gives costs ≤
$$
3(y+2)(3y^2/8+4y+15)e^{24-11y/16}.
$$
For example the pointwise second derivative of Ze^{cp}T_y is bounded by
$3(y+2)(c^2+6c+y+24)e^{-11y/16}$; the first derivative costs at most the same factor without the bracket times c+3. Add this to one half the c² and second derivative costs and use c≤y/2. Both displayed polynomial times exponential bounds decrease on y≥80: each polynomial has nonnegative coefficients and degree at most three. Divide by c x ≥28·16 by using g₀≥xe^{-x}. Substituting y=80 gives
$$
\frac{|E'(0)|+\int_0^{1/2}|c^2E-E''|}{c g_0}<1/1000.
$$

For the high spectrum use 0≤U,V≤e^{-pm}/p and the convexity and monotonicity already stated. In particular ∫_{1/2}^{1-κ}|D|≤2e^{-m/2}/m, and
$$
\int_{1/2}^{1-κ}|D''|+|D'(1-κ-)|\le-U'(1/2)-V'(1/2)\le(4m+24)e^{-m/2}.
$$
For the first inequality apply nonnegative U'', V'' and the fact the terminal one-sided slopes are nonpositive; at that endpoint |U'−V'|≤−U'−V'. Thus the high-spectrum cost divided by cg₀ is at most
$$
\left[2c/m+(4m+24)/c\right] e^{(x-c)/2}/x<1/10,
$$
since m/c≤13/7 and x−c≤−4. The constant tests, including the elementary e bounds, can be independently re-run with [exact rational arithmetic](../../certificates/verify_outer_moderate_constants.py).

### 2.3 Signed transform

Extend G affinely past 1/2. The cosine transform of its leading constant times e^{-cp} is g₀c/(c²+ω²). The linear first-slope term, after multiplication by c²+ω², costs at most |G'(0)|. For the curvature term use the convolution with s e^{-cs}; its complex Fourier multiplier is (c+iω)^{-2}, so the cost is at most ∫₀^{1/2} e^{-cp}|G''| ≤192g₀/c. Subtract the affine extension from the physical D. The residual is E below 1/2 and D minus the affine extension above. There is no value or first derivative jump at 1/2, and the only physical cutoff jump is the true first-derivative atom. Two integrations by parts with c²−∂² bound the other costs by the low and high quantities in §2.2 and the affine-tail cost 2|G'(1/2)|e^{-c/2}. Dividing by cg₀ gives the strict bound
$$
1-\frac{1}{112}-\frac{192}{28^2}-\frac1{1000}-\frac1{10}-\frac1{1000}>1/2.
$$
All real frequencies therefore have positive physical cosine transform for this unbounded branch. Combined with §1 and the reconstructed moderate bridge, this completes 16≤x≤24 for all y≥32 when x≤N. No claim concerning x<16 or the remaining deep near-diagonal branch follows from this document.
