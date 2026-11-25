# The full continuum deep near-diagonal strip

**Theorem.** Put
$$
 J_v(p)=\int_v^\infty \frac{\sinh((1-p)z)-(1-p)z}{\cosh z-1}\,dz,
 \qquad \Delta=J_m-\sqrt{J_{m-c}J_{m+c}}.
$$
Use the removable regularization at $p=0$ and extend by zero past $p=1$. If $x=m-c\ge24$ and $0<c\le9$, then the half-line cosine transform of $\Delta$ is **strictly positive at every real frequency**. All estimates are uniform down to $c=0$ after division by $c^2$; at $c=0$ the unnormalized defect is zero. This is rigorous interval-assisted mathematics for **$\kappa=0$**, not a finite-$\kappa$ theorem. Together with [the separation proof](deep-separation-proof.md), which admits $\kappa=0,\ c\ge9$, this exhausts the continuum deep region without an ineffective radial threshold.

## 1. Exact splitting

Put $u=1/m,\ h=c/m,\ q=1-p,\ T=mq$. Then $0<u\le1/(24+c)$ and $0\le h\le3/11$. Define
$$
 M=q[\cosh(cq)-1]/c^2,\quad N=q\sinh(cq)/c,\quad
 a=(\cosh c-1)/c^2,\quad A=(1+u)M-uN,\quad L=u e^{-mp}/(1+p).
$$
Extend $M,N$ by zero after $1$ and keep $L$ on the entire half-line. Let $W_v=pe^{pv}J_v$ and $g_v=\sqrt{W_v}$. Identity (2) of [the expansion](deep-finite-kappa-asymptotic.md) gives on $0<p<1$:
$$
 \frac{\Delta}{2m e^{-m}c^2}=A-aL+u^3(Z+B)+R,\qquad
 B=\frac{e^{-T}}{1+q}\left(\frac{\sinh(hT)}h\right)^2,\qquad
 Z=p\left\langle e^{T/2}\frac{\partial_\rho g_{m\rho}}p\right\rangle^2. \tag{1}
$$
Here $\rho=1+sh$ and $\langle k\rangle=\frac12\int_{-1}^1 k(1+sh)\,ds$. For the square use $(g_{m+c}-g_{m-c})/(2h)=\langle\partial_\rho g_{m\rho}\rangle$, so no interval difference is divided by a vanishing half-gap. Explicitly, the three linear columns at alias $j\ge1$ are
$$
 e^{-(j-1)m}q[(1+u/j)G(j-p)-uH(j-p)],\quad
 \frac{ju e^{-(j-1+q)m}}{2(j+q)}G(j-1+2q),\quad
 -\frac{(j+1)u e^{-(j-1+p)m}}{2(j+p)}G(j),
$$
where $G(t)=(\cosh(ct)-1)/c^2,\ H(t)=\sinh(ct)/c$. For $j=1$ these are $A,u^3B,-aL$; only $j\ge2$ remain in $R$. The narrow negative spike is **not** included in an absolute curvature estimate.

On $[0,1]$, each of $|R|,|R'|,|R''|$ is less than $1/250$. In fact $G(t)\le t^2e^{ct}/2,\ H(t)\le te^{ct},\ G'=H,\ G''=\cosh(ct)$ for $t\ge0$. Either derivative through order two of the first column costs at most $4(j+1)^2 e^{-(j-1)m+cj}$, also using $|H''|\le9e^{ct}$. In the second and third columns the rational $j/(j+q),(j+1)/(j+p)$ and their first two derivatives cost at most 2, and $G$ with its affine argument and its first two derivatives cost at most $(j+1)^2$ times the exponential of $c$ times that argument. By the product rule and $u(m+1)^2\le2m$ each costs at most $20m(j+1)^2$ times its corresponding mixed exponential. Its exponent in the second is $-(j-1)(m-c)-q(m-2c)\le-(j-1)(m-c)$, and in the third at most $-(j-1)(m-c)+c$. Since $m e^{-(j-1)(m-c)}\le33e^{-24(j-1)}$, the summed bound is
$$
 1324\sum_{j\ge2}(j+1)^2 e^{-24(j-1)+9}<1/250.
$$
This also justifies termwise differentiation. Beyond $p=1$ prescribe $Z=B=0,\ R=aL$. At cutoff the true defect is $q^3$ times a smooth function for fixed radii: extract $q^3$ in the ordinate integrals, leaving smooth strictly positive factors. The model terms also vanish there with zero first derivative (for $B,Z$, see below). There is thus no artificial slope jump of $R$. Past 1, $\int a L<e^{-24}$ and $\int aL''=-aL'(1)<2e^{-24}$, since $a\le e^c/2$. Consequently $R$ costs at most $41/10000$ in the ordinary transform, or $81/10000$ in the $\omega^2$-scaled transform including the initial slope.

## 2. Two nonsingular charts for the square

In this section derivatives are in $T$ with $m,c,\rho$ fixed. Let $v=m\rho\ge24,\ t=\rho T,\ v_0=1/v,\ z=e^{-v},\ \ell=j+1,\ b_j=v_0/\ell$. Distinguish $t$ and $T$. The usual expansion $(\cosh w-1)^{-1}=2\sum_{\ell\ge1}\ell e^{-\ell w}$ gives
$$
 W=\sum_{j\ge0}z^j\left[\frac{p\ell}{j+p}-
 \frac{p\ell e^{-2t}}{\ell+q}-2p(t+q/\ell)e^{-t}\right],\qquad
 H_v=e^{pv}F_p(v)=\sum_{j\ge0}\ell z^j(1-e^{-2t}-2t e^{-t}).
$$
The first summand inside the bracket for $j=0$ is 1. Differentiation at the lower limit gives the exact radial ODE $\partial_\rho g_{m\rho}/p=m(W-H_v)/(2\sqrt W)$. Put $E(t)=(\sinh t-t)/t^3,\ C(t)=(\cosh t-1)/t^2$ with removable values.

**$T\le2$.** Direct rearrangement at each alias gives
$$
 W=t^3K,\qquad K=2p e^{-t}\sum_{j\ge0}z^j
 \frac{E+b_j C+b_j^2+b_j^3}{1-(q/\ell)^2}.
$$
Cancel $p$ in the first term to obtain
$$
 K=\frac{2 e^{-t}(E+v_0 C+v_0^2+v_0^3)}{1+q}+\epsilon_0,\qquad
 e^{T/2}\frac{\partial_\rho g}p=T^{3/2} k,
$$
$$
 k=\frac{e^{T/2-t}\rho^{3/2}}{\sqrt K}
 \left[\frac{C/\rho+u/\rho^2+u^2/\rho^3-TE}{1+q}+\epsilon_1\right]. \tag{2}
$$
Thus $Z=pT^3\langle k\rangle^2$ is explicitly smooth at $T=0$. Each Taylor lane (function, derivative, half second derivative) of $\epsilon_0,\epsilon_1$ is enclosed in $\pm10^{-7},\pm10^{-6}$ respectively.

Indeed use the submultiplicative jet norm $\|f\|_J=|f|+|f'|+|f''|/2$. Here $q\le1/12,\ \rho\le14/11<1.3,\ t\le2.6,\ \cosh(2.6)<7$. The integrals $C=\int_0^1(1-s)\cosh(st)ds,\ E=\frac12\int_0^1(1-s)^2\cosh(st)ds$ give $\|C\|_J<5.6,\|E\|_J<1.7$. Also $\|(1-q^2/\ell^2)^{-1}\|_J<1.1,\ b_j\le1/24$. Put $A_j=(E+b_j C+b_j^2+b_j^3)/(1-q^2/\ell^2)$, so $\|A_j\|_J<2.2$. The two errors are $2p e^{-t}\sum_{j\ge1}z^j A_j$ and $m\sum_{j\ge1}z^j(pA_j-\ell E)$. Using $m z\le33e^{-24}$, the rational geometric sums give $16 e^{-24}$ and $7(33)e^{-24}$, respectively.

**$T\ge2$.** With $D=e^{-t}$, the same formulas, without subtracting two near-unit values, are
$$
 W=1-\frac{pD^2}{1+q}-2p(t+q)D+\epsilon_2,
$$
$$
 e^{T/2}\frac{\partial_\rho g}p=
 \frac{e^{T/2-t}T[D/(1+q)+t-p]
 +b_*[-1/(1+p)+(1+2q)D^2/(2+q)+((1+q)t-pq/2)D]+\epsilon_3}{\sqrt W}. \tag{3}
$$
Here $b_*=m e^{-v+T/2}$. The term with this factor is the exact $j=1$ contribution to $m e^{T/2}(W-H_v)/2$. With $M_0=\max(24,\inf T)$ on any closed $T$-cell, $b_*\le\min(33e^{-24+T/2},M_0 e^{-M_0+9+T/2})$, since $v\ge24,\ v\ge m-9,\ m\ge T$. Its derivative lanes are $b_*,b_*/2,b_*/8$.

Each lane of $\epsilon_2,\epsilon_3$ costs at most $10^{-7}$ and $100 b_*e^{-24}$. Indeed $\|p\|_J,\|q\|_J<1.05,\ \|D\|_J<3.15,\ \|D^2\|_J<7.1,\ \|tD\|_J<5.3$. The jet norms of the reciprocals $1/(j+p),1/(\ell+q)$ are at most $1.05/j,1.05/\ell$. For $j\ge1$, the three terms of $W_j/z^j$ cost less than 25, namely at most $2.21+7.83+2.1(5.3+1.05\cdot3.15/2)<25$. Summing gives $26e^{-24}$. On removing $z^j$, the summands of $(W-H_v)/2$ are
$$
 \frac{\ell}{2}\left[-\frac j{j+p}+\frac{j+2q}{\ell+q}D^2+
 2\left(\frac{j+q}{\ell}t-\frac{pq}{\ell^2}\right)D\right].
$$
For $j\ge2$, their norm is at most $13\ell$: the four terms inside the brackets cost respectively $1.05,\ 1.5(7.1),\ 2(1.03)(5.3),\ 2(1.05)^2(3.15)/9$ or less. Multiplying by $e^{T/2}$ costs $13/8$ times its value. Since $\sum_{j\ge2}\ell z^j<3.01z^2$, the claimed $100 b_* e^{-24}$ bound follows.

An independent derivation of every input-alias allowance and the unbounded tail estimate is in the endpoint audit.

## 3. Compact enclosure of an unbounded-depth layer

Also $B=\langle T e^{(sh-1/2)T}\rangle^2/(1+q)$, valid even at $h=0$. On $0\le T\le96$, [the layer verifier](../../certificates/deep_continuum_layer.py) uses (2) and (3) to bound $|Z|+|B|$ and $|Z''|+|B''|$. The $s$-integrals use four equal closed cells with exact weights. Two closed cells cover $u$ in each of nine unit $c$-bands. Intersect $q=uT$ with $[0,1]$, $u$ with $u\le1/\inf T$ when $\inf T>0$, and $h$ with $h\le c_1/(24+c_1)$ for the upper band endpoint $c_1$. This does not differentiate a clamped function: on the true parameter pieces $q'=u$, as used in the derivative lanes. Empty cells are omitted. Every Taylor lane is evaluated over the full spectral interval, not only at its midpoint. The series of $E,C$ and their two derivatives have positive terms, and tails after index 16 less than $10^{-22}$ for $t\le3, |t'|\le1.3$ (the successive ratios are less than $1/2$).

For $T\ge96$, use $m\ge T,\ |t-T|\le9,\ v\ge m-9,\ q\le1,\ u\le1/96,\ \rho\le1.1$. Each of $W-1,W',W''/2$ is at most $10^{-5}$ in absolute value by (3), giving $\|1/\sqrt W\|_J<2$. The first numerator of (3) costs at most $16(T+12)^2 e^{9-T/2}$ since $\|e^{T/2-t}\|_J<2e^{9-T/2},\ \|T\|_J=T+1$ and $\|D/(1+q)+t-p\|_J<T+12$. The second and error together cost at most $40T e^{9-T/2}$: $b_*\le T e^{9-T/2}$ since $m\ge T$, and the explicit $j=1$ bracket in (3) has norm less than 23 (here $t\ge87$). Thus the jet norm of the averaged ordinate in (3) is $\le1000(T+12)^2 e^{9-T/2}$; the one in $B$ costs $\le5(T+1)e^{9-T/2}$. Product differentiation bounds both the summed masses and summed absolute second derivatives by $2\cdot10^7(T+12)^4e^{18-T}$. Its integral above 96 is $<10^{-8}$ by four integrations by parts.

Now $dp=-u\,dT$ and each $p$ derivative multiplies a $T$ derivative by $-m$. Thus $u^3(Z+B)$ has mass at most $u^4\int(|Z|+|B|)$ and total absolute second derivative at most $u^2\int(|Z''|+|B''|)$ (integrals in $T$). At $p=1$ both the values and the first derivatives vanish by (2); hence the initial slope at zero is bounded in absolute value by the curvature budget as well. For the bandwise upper bounds $M_b,C_b$ this costs $\min(\omega^2 M_b,2 C_b)$ in the $\omega^2$-scaled transform. The larger interval overestimate near $c=9$ is harmless because the model there has a larger surplus. [The table](../../certificates/deep_continuum_tables.py) is **verified** by the layer computation before being used in the frequency proof.

## 4. All frequencies of the main term

Use $\widehat f(\omega)=\int_0^\infty f(p)\cos(\omega p)\,dp$. The positive mixture $\widehat L=u\int_0^\infty e^{-s}(m+s)/((m+s)^2+\omega^2)\,ds$ gives $\widehat L\le u^2$ and $\omega^2\widehat L\le(1+u)\omega^2/(m^2+\omega^2)$. Since $A,A'$ vanish at 1,
$$
 \omega^2\widehat A=\int_0^1 A''(p)(1-\cos(\omega p))\,dp,\qquad
 A''=(1+u)q[2\operatorname{sinhc}(cq)+\cosh(cq)]-u[2\cosh(cq)+(cq)^2\operatorname{sinhc}(cq)].
$$
There is no assumed sign of $A''$. For $0\le|\omega|\le1$, divide by $\omega^2$ and enclose the average test weight on a cell between $\overline{p^2}(1/2-p_{\rm hi}^2/24)$ and $\overline{p^2}/2$, also valid by a limit at zero. For $1\le|\omega|\le32$, use closed frequency intervals of width $1/4$ and $p$ intervals of width $1/64$. At center frequency $w$ the true cell average of the test is $1-\cos(w p_{\rm mid})\operatorname{sinc}(w/128)$. Its variation is at most $p_{\rm mid}/8$ since $|\partial_w\cos(wp)|\le p$. Trigonometric enclosures use bisection to a point argument at most 1, alternating series with the next term included, then the double-angle formulas for cosine and sinc. The width of the $c$ intervals is $1/8$.

For $|\omega|\ge32$, $M'',N''$ are both nonnegative and nonincreasing in $p$. For each, integration by parts bounds the absolute integral against cosine by its value at $p=0$ divided by $|\omega|$, including any nonzero boundary value at $p=1$. After using $\omega^2a\widehat L\le(1+u)a$ and the initial slopes of $M,N$, a valid lower bound for $\omega^2(\widehat A-a\widehat L)$ is
$$
 \operatorname{sinhc}c-u\cosh c-
 \frac{(1+u)(2\operatorname{sinhc}c+\cosh c)+u(2\cosh c+c^2\operatorname{sinhc}c)}{32}.
$$
See also an independent short high-frequency audit. Deduct the branch-specific costs of §§1–3. [The frequency verifier](../../certificates/deep_continuum_frequency.py) checks all nine groups. The least remainders in its units on $[0,1],[1,32],[32,\infty)$ respectively exceed $97522,\ 67606,\ 715719$ divided by $10^6$ (the last two are $\omega^2$-scaled). It uses integer endpoints on the $2^{-112}$ grid, rigorous integer square root, and range-reduced exponential.

**Complete replay:** run the following command from the project root with Python 3, normal (non-optimized) mode.

    sh certificates/replay_deep_continuum.sh

It checks the numerical majorants, rederives every table entry, then verifies all model signs and deductions. [Complete output：从源码重放](../../replay.sh). This certificate closes the all-frequency continuum strip and is not a finite-$\kappa$ interval proof or an unconditional Lean theorem.

