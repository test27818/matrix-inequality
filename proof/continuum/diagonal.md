# 同半径余量与零频谱端点

## 3. Continuous Fourier transform and the open mixed sign

For $A,B>0$, set $K=(A+B)/2$ and

$$b_{A,B}(t)=e^{-Kt}/t-\mathbf1_{0<t<1}\sqrt{J_{1-t}(A)J_{1-t}(B)},\qquad D(A,B,\omega)=\int_0^\infty b_{A,B}(t)\cos(\omega t)\,dt.$$

On the diagonal $A=B$, direct differentiation under the integral and the elementary cosine transform below yield an explicit, nonnegative slack:

$$D(A,A,\omega)=\frac{1-\cos\omega}{\omega^2}\int_A^\infty\frac{2v^3e^{-v}}{(1-e^{-v})^2(v^2+\omega^2)}\,dv,$$

where the factor at $\omega=0$ is $1/2$. To verify it, integrate $[\sinh(A(1-t))-(1-t)A]\cos(\omega t)$ from $0$ to $1$: the sine-hyperbolic part integrates to $A(\cosh A-\cos\omega)/(A^2+\omega^2)$, and the linear part to $A(1-\cos\omega)/\omega^2$. The free exponential derivative contributes $-A/(A^2+\omega^2)$. Consequently $-\partial_A D=(1-\cos\omega)A^3/[\omega^2(A^2+\omega^2)(\cosh A-1)]$. The integration constant vanishes at $A=+\infty$, as is also seen from the finite-$n$ nonnegative slack formula above, or by the cancellation at $t=0$ in the defining transform.

For unequal radii the exact leftover is

$$D(A,B,\omega)-D(K,K,\omega)=\int_0^1\left[J_{1-t}(K)-\sqrt{J_{1-t}(A)J_{1-t}(B)}\right]\cos(\omega t)\,dt.$$

Thus the diagonal term is rigorously positive away from nonzero $2\pi\mathbb Z$; at those phases all mixed positivity must come from this last transform. Pointwise sign of its amplitude alone would not imply a nonnegative Fourier transform at arbitrary frequencies. This document does not claim the general mixed sign or the unrestricted all-orders formalization.


### Large frequencies for each radial pair, uniformly on compact radial ranges

Write $c(s)=s\coth(s/2)$ and $b_s=b_{s,s}$. For $0<p<1$, an integral with no singularity at $p=0$ is
$$
 b_s(p)=\int_s^\infty
 \frac{(1-p)v+e^{-v}\cosh(pv)-e^{-pv}}{\cosh v-1}\,dv.
$$
The integrand and every order of its $p$ derivatives are integrably dominated for, say, $|p|<1/2$ and $s>0$. Thus $b_s(0)=2s/(e^s-1)$ (integrate, or compare $s$ derivatives and the limit at infinity) and $b_s'(0)=0$ (the derivative of the integrand vanishes pointwise). Consequently
$$
J_{1-p}(s)=p^{-1}-c(s)+(s^2/2)p+O(p^2).
$$
This also justifies the zero integration constant for the diagonal transform: the identity for $b_s$ above gives $\int_0^1|b_s(p)|dp\to0$ as $s\to\infty$ by domination, and $\int_1^\infty e^{-sp}/p\,dp\to0$.

Define $\Delta_{A,B}(p)=J_{1-p}(K)-\sqrt{J_{1-p}(A)J_{1-p}(B)}$ on $0<p<1$, and let it be zero for $p>1$. Taylor expansion gives
$$
-\Delta'_{A,B}(0)=\frac{(A-B)^2-(c(A)-c(B))^2}{8}=:\delta_{A,B}.
$$
Here $0<c'(s)<1$ for $s>0$: $c'(s)>0$ follows from $\sinh s>s$, and $c'(s)<1$ from $1-e^{-s}<s$. Near $p=1$, the convergent Taylor series for $J_\mu$ factors as $\mu^3$ times a strictly positive analytic function of $\mu^2$. Hence $\Delta,\Delta',\Delta''$ vanish at 1, and $\Delta'''$ is integrable. Twice integrating by parts yields
$$
\omega^2 D(A,B,\omega)=\omega^2 D(K,K,\omega)+\delta_{A,B}-\int_0^1 \Delta_{A,B}''(p)\cos(\omega p)\,dp.
$$
Once more by parts bounds the final integral in absolute value by $\|\Delta'''\|_{L^1}/|\omega|$. Thus for every unequal pair $D>0$ at all $|\omega|>\|\Delta'''\|_{L^1}/\delta_{A,B}$, and $\omega^2D\to\delta_{A,B}>0$. The diagonal term above is explicitly nonnegative and $O(\omega^{-4})$.

