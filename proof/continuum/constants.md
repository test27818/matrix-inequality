# Direct scalar constant identity, with no finite-order Hermite input

For the continuous coefficients in [the matrix transfer](sampling.md), individual roots-of-unity contacts are not needed for the matrix sum. The following elementary identity suffices for every $n\ge2$ and $x>0$:
$$
 (n-1)g_0(r)-\sum_{j=1}^{n-1}g_j(r)=-\log S_n(r),
 \qquad r=e^{-x/n},\quad S_n(r)=1+r+\cdots+r^{n-1}.
$$
Here $g_0(r)=-x/[n(e^x-1)]$ and, for $1\le j<n$,
$$
 g_j(r)=\frac1n\int_x^\infty
 \frac{\sinh((1-j/n)v)-(1-j/n)v}{\cosh v-1}\,dv .
$$
Indeed the transfer defines $g_0=-T_x(0)/(2n)$. At $\kappa=0$ the integrand of $T_x(0)$ equals
$2e^{-v}(v-1+e^{-v})/(1-e^{-v})^2=-\partial_v[2v/(e^v-1)]$,
so that $T_x(0)=2x/(e^x-1)$ exactly.

**Proof.** Put $H(x)=-(n-1)g_0(r)+\sum_{j=1}^{n-1}g_j(r)$. Differentiate the displayed coefficients, substituting $j\mapsto n-j$ inside the finite sum. Since $\sum_{j=1}^{n-1}j/n=(n-1)/2$ and
$$
 \frac1{2(\cosh x-1)}=\frac{e^x}{(e^x-1)^2},
$$
the terms proportional to $x$ cancel. Therefore
$$
 n H'(x)=\frac{n-1}{e^x-1}
       -\frac{\sum_{j=1}^{n-1}\sinh(jx/n)}{\cosh x-1}.
$$
The quotient on the right equals $1/(e^{x/n}-1)-1/(e^x-1)$. For a check, put $t=e^{x/n}>1$; summing the ordinary geometric progressions for $t^j$ and $t^{-j}$ shows that the quotient is $(t^n-t)/[(t^n-1)(t-1)]$. Consequently
$$
 H'(x)=\frac1{e^x-1}-\frac{1}{n(e^{x/n}-1)}
 =\frac{d}{dx}\log\frac{1-e^{-x}}{1-e^{-x/n}}.
$$
Both sides of the undifferentiated identity tend to zero at infinity: each $J_x(j/n)$ is a tail integral of an absolutely integrable function and $x/(e^x-1)\to0$. This proves the identity.

After the scalar pairwise lower bound has been established by sampling (or directly), its constants sum to $(n-1)\sum_i g_0(r_i)$, and positive semidefiniteness of the Schur powers bounds all off-diagonal nonconstant modes below by $-\sum_{i,j\ge1}g_j(r_i)$. Substituting the identity above completes the matrix assembly. Thus all-order root lifting, convergence of the finite-order coefficients, and individual contacts can be removed from this route. This shortcut concerns transfer only and does not by itself establish mixed scalar positivity on a new region.
