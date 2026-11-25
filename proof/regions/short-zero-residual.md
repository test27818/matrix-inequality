# A bounded shallow-to-medium rectangle for the original spectrum

For $0\le x\le3$ we need the midpoint gap $D=J_m-\sqrt{J_xJ_y}$, $m=(x+y)/2$, to have a nonnegative cosine transform for every frequency. This is the *finite parameter* spectrum with $b=\kappa^2\in[0,1/256]$. One takes $p\in[0,1-\kappa]$, $q=1-p$ and extends $D$ by zero on the right. We divide the spectrum at $p=1/2$; there is no artificial atom there. The key is that even at $t=0$ the singular term at $p=0$ occurs in just one exponential alias. Removing it gives two convergent positive coefficient evaluations and no radial quadrature.

## Spectral expansion and bounds

Write $P=q(q^2-b)$ and $S_n=\sum_{h=0}^n q^{2h} b^{n-h}$. Let
$$
 C_{j,n}(t)=\frac{2e^{-jt}}{j^{2n+3}}\sum_{r=0}^{2n+3}\frac{(jt)^r}{r!},
 \qquad R_t=\sum_{n\ge0}S_n\sum_{j\ge2}C_{j,n}(t),
 \quad A_t=R_t+\sum_{n\ge0} S_n C_{1,n}(t).
$$
These come from expanding $1/(\cosh v-1)=2\sum_{j\ge1}j e^{-jv}$ and $(\sinh(qv)-q\sinh(\kappa v)/\kappa)/P=\sum_{n\ge0}S_n v^{2n+3}/(2n+3)!$; the second identity extends to $q=\kappa$, and to $\kappa=0$ using the continuous ratio. On the physical open interval $J_t=P A_t$. Each $C_{j,n}$ is nonincreasing in $t$, since it is $2j$ times the corresponding positive exponential moment on $[t,\infty)$.

**Low part.** Set $H_t=2\sum_{n\ge0}b^n e^{-t}\sum_{r=0}^{2n+1}t^r/r!$. Integrating the whole $j=1$ alias exactly gives
$$
 T_t:=e^{-pt}/p-J_t =\frac{e^{-(1+q)t}}{1+q}+ qH_t-P R_t .
$$
The right-hand expression is smooth also at $p=0$. Define $W_t=1-pe^{pt}T_t$ and $c=(y-x)/2$; factor the difference to remove $1/p$:
$$
 D=e^{-cp}\frac{T_x}{1+\sqrt{W_x}}-T_m+
  e^{cp}\frac{\sqrt{W_x}T_y}{1+\sqrt{W_y}} .
$$
For $p\le1/2$ the verifier uses this expression and its first two derivatives. The derivative variables are always $p$, including when evaluating series in $q=1-p$. Positive lower bounds on square root denominators are tested in each cell.

**High part.** For $q\le1/2$ use $D=P(A_m-\sqrt{A_xA_y})$ and two derivatives from the positive representation. A mesh cell straddling $q=\kappa$ also encloses zero in its range for $D,D''$, because of the extension. The boundary mass of the distribution $D''$ at $p=L=1-\kappa$ is $-D'(L-)=2b(A_m-\sqrt{A_xA_y})$. We include its independently computed lower bound capped at zero; at $q^2=b$ the exact polynomial is $S_n=(n+1)b^n$. Hence the proof does not assume this atom to have the favorable sign.

**Two truncations.** With $N=36$, all jets retain coefficients $0\le n\le N$. For $q\le1/2$, $S_n,-\partial_pS_n,\partial_p^2S_n/2$ are nonnegative and bounded by $16(n+1)^3 4^{-n}$: for each of the $n+1$ terms one can replace $q$ by $1/2$, $b$ by $1/4$, losing at most a factor $4(2n)^2$ in two differentiations. Also $\sum_{j\ge1}C_{j,n}\le 2\zeta(3)<3$. For $R$ alone at $q\le1$,
$$
 \sum_{j\ge2}C_{j,n}\le2\left(2^{-2n-3}+\int_2^\infty s^{-2n-3}\,ds\right)\le \tfrac12\,4^{-n},
$$
and the three polynomial lanes are bounded by $4(n+1)^3$. The tail geometric ratio of $(n+1)^3/4^n$ is below $1/2$ for $n\ge37$. Thus each omitted lane in both regimes has absolute error at most $96(N+2)^3 /4^{N+1}< 2^{-41}$, with rigorously known lane signs. The remainder for $H$ before the factor 2 is at most $b^{N+1}/(1-b)$, as a Poisson partial exponential sum is at most 1.

For each $n$ and a particular endpoint $t$, summing only $1\le j\le M$ costs at most $1/((n+1)M^{2n+2})$ by the decreasing integral test. If $t>0$ one can also use $C_{j,n}\le j e^{-(j-1)t}C_{1,n}$ from the positive integral representation and take the smaller of that previous upper bound and $C_{1,n}(M+1)e^{-Mt}/(1-e^{-t})^2$. The program selects any convenient $M\ge3$. The monotonicity in depth encloses each coefficient throughout the closed parameter box by its two endpoint evaluations.

## Why positive budgets imply all real frequencies

We enclose $D,D''$ on every closed spectral interval and compute a rigorous lower bound $a\le\min(0,-D'(L-))$. The result follows either from the [signed-curvature frequency bounds](signed-curvature.md) or from the [Helmholtz bound](moderate-helmholtz.md): a multiplier $\lambda>0$ suffices when
$$
 \int_0^1\min(\lambda^2 D,2D''-\lambda^2D)\,dp+2a>0 .
$$
Here the bulk derivative in the integrand excludes its cutoff atom. The nine or fewer choices of $\lambda$ are *alternative* sufficient bounds, not assumptions on the frequency. More generally the elementary two integrations by parts give $\omega^2\int D\cos=\int D''(1-\cos)+(-D'(L-))(1-\cos(\omega L))$; using the piecewise lower bounds for the nonnegative weights and their high-frequency variation bound yields the alternate signed-frequency test. These tests, the dyadic depth and phase subdivisions, and the include-zero treatment cover continuous parameters, not only sample points.

## Reproducibility

The verification source is [short_zero_residual.cpp](../../certificates/short_zero_residual.cpp). The enclosure endpoints are binary64 dyadics; every arithmetic endpoint is extended to the next representable value in the outward direction. The host must use IEEE-754 binary64 round-to-nearest, gradual underflow, and no unsafe floating optimizations. No transcendental library values decide an interval: exponentials use range reduction, 18 Taylor terms on $[0,1/2]$ and a proved geometric remainder, and sine/cosine weights use alternating series and doubling; the square root only proposes endpoints and tests their squares with outward operations. Unsatisfied preconditions, overflow, and failed subdivisions cannot count as a successful proof. Formula cutoffs, proof status, and root exact positions use dyadic integers. Underflow to zero is enclosed by the adjacent subnormal; compilation should use the flags below.

    c++ -O3 -std=c++20 -ffp-contract=off -fno-fast-math -fsanitize=undefined \
        -fno-sanitize-recover=all certificates/short_zero_residual.cpp -o /tmp/szero
    /tmp/szero 0 3 8 24 32 19 40000
    /tmp/szero 0 3 6 8 64 22 80000
    /tmp/szero 0 1 4 24 32 19 100000

Certified counts are in [short-zero-residual-replay.txt：从源码重放](../../replay.sh); an unfinished outer box never counts. For the argument from the continuum of real frequencies to the finite Hermite kernel see [the Fourier sampling lemma](../continuum/positive-periodization.md).

**Final complete rectangle.** The integer signed-curvature checker also passes the whole triangle $0\le x\le3$, $x\le y\le6$ in one independent run with *its own conservative atom calculation*, not even requiring a separate positivity lemma:

    PYTHONDONTWRITEBYTECODE=1 python3 certificates/run_signed_curvature.py 6 3 8 16 6

[Full integer replay：从源码重放](../../replay.sh): the last level has 23,924 pending radial pairs; it finishes with zero unresolved. Combining this with the residual runs $6\le y\le8$ (682 accepted leaves) and $8\le y\le24$ (508) proves every pair in $0\le x\le3$, $x\le y\le24$, for the original finite $0\le\kappa\le1/16$ and every real frequency.

The two nonredundant residual runs were also repeated with optional `PSD_ZERO_TRACE` filenames [short-zero-6-8-cells.tsv：从源码重放](../../replay.sh) and [short-zero-8-24-cells.tsv：从源码重放](../../replay.sh). Each entry gives the accepted closed three-parameter box, spectral mesh, multiplier and the (dyadic) strict lower margin. All 1,190 of these leaves used the elementary Helmholtz option, so no trigonometric certification was needed in these two sets. The independent exact geometry check verifies interior disjointness and the total volume covering both root boxes; [its output：从源码重放](../../replay.sh) also lists the smallest certified lower margin of each file. Leaf verification remains the job of the source verifier.
