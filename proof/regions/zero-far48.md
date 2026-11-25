# Complete zero-depth far strip, finite kernel: 0 ≤ x ≤ 3, y ≥ 48

**Theorem.** For the original physical $J_v$ and
$D=J_m-\sqrt{J_xJ_y}$, with $0\le\kappa\le1/16$, $L=1-\kappa$,
$m=(y+x)/2$, $c=(y-x)/2$, one has for all real $\omega$
$$
 (c^2+\omega^2)\int_0^L D(p)\cos(\omega p)\,dp
 \ \ge\ \frac{35957}{81000}c g_0\ >\ \frac25 c g_0 ,
 \qquad g_0=T_x(0)/2>0 .
$$
Here $T_t=e^{-pt}/p-J_t$ and $W_t=1-p e^{pt}T_t$ have their regularized endpoint values. We use no limiting substitution for $\kappa=1/N$.

## 1. A small closed rational jet check

Put $b=\kappa^2$, $q=1-p$, $z=q^2$, $P_0=q(z-b)$ and $C(v)=(\cosh v-1)/v^2$ with $C(0)=1/2$. Use the removable positive series from [the shallow square](finite-shallow-square.md),
$$
 f(q,v)=\frac v{C(v)}\sum_{a,h\ge0}
 \frac{(qv)^{2a}(\kappa v)^{2h}}{(2a+2h+3)!},\quad U_x=\int_0^x f(q,v)\,dv ,\quad F_p=P_0 f .
$$
Integrating $1/(\cosh v-1)=2\sum_{j\ge1} j e^{-jv}$ gives
$\int_0^\infty f=2/[(1-z)(1-b)]+C_*(q)$ where
$$
 C_*(q)=2\sum_{j\ge2}\frac{j}{(j^2-z)(j^2-b)} .
$$
This holds on $p>0$; cancelling the $j=1$ pole explicitly gives a regular formula valid also at $p=0$:
$$
 T_x=\frac{2z+2q+1-b}{(1+q)(1-b)}
       -P_0 C_*(q)+P_0 U_x-\int_0^x e^{-pv}\,dv .             \tag{1}
$$
In particular this identity has no $1/p$ or $1/x$ and works at $x=0$. The tail after $j=n=48$ of $C_*$ and its first and second $p$ derivatives lies in
$[0,M],[-2M/u,0],[0,(2/u+8/u^2)M]$, respectively,
where $u=(n+1)^2-1$ and $M=((n+1)^2/((n+1)^2-b))/(n(n+1))$.
Indeed $2j/[(j^2-z)(j^2-b)]\le2((n+1)^2/((n+1)^2-b))/[j(j^2-1)]$
and $\sum_{j>n}1/[j(j^2-1)]=1/[2n(n+1)]$; differentiating the other denominator yields the two displayed factors with these signs.

Let $G=T_x/(1+\sqrt{W_x})$. On the full box $(x,p,b)\in[0,3]\times[0,1/2]\times[0,1/256]$,
[the exact checker](../../certificates/verify_zero_far48.py) gives
$$
 W_x\ge9/64,\quad g_0=G(0)\ge1/16,\quad
 |G'(0)|\le3/5,\quad |G'(1/2)|\le1,\quad |G''|\le4,\quad
 |(\sqrt{W_x})'|\le3,\quad |(\sqrt{W_x})''|\le21 .           \tag{2}
$$
Unnormalized derivatives are used in every jet. For (1) the program integrates positive radial cells at spacing $1/32$ for $f,f_q,f_{qq}$ and $e^{-pv},v e^{-pv},v^2 e^{-pv}$; the two first $p$ derivatives of the integrated $q$-profile have signs $-,+$, respectively. On the cell ending at $x_{i+1}$ and beginning at $x_i$, the lower integral through $x_i$ and upper integral through $x_{i+1}$ bound every intermediate endpoint. From the defining positive series, the $j$-th $q$ derivative of $f$ ($j=0,1,2$) lies between $v^{j+1}E^{(j)}(qv)/C(v)$ and that expression times $(1-bv^2/20)^{-1}$, where $E(t)=\sum_{a\ge0}t^{2a}/(2a+3)!$; compare termwise, including after differentiation, and use $(2a+3)!/(2a+2h+3)!\le20^{-h}$. All $b$ are enclosed together with $b$ replaced by the interval $[0,1/256]$. This also encloses the lower factor $1$ exactly or outward. The positive coefficient Taylor enclosures for $E,E',E'',C$ use 48th degree and the verified tail $2\cdot4^{49}/49!<2^{-56}$. For $e^{pv}$ the Taylor polynomial through degree 24 is bounded above by adding $2\cdot2^{25}/25!<2^{-56}$ (since $pv\le3/2$). Reciprocals and square roots are enclosed by integer division and integer square root. There is no physical square-root clamp: the lower bound on $W_x$ is checked **before** the square root. The first and last of 32 closed $p$-bins also cover the endpoint slope bounds. The [replay：从源码重放](../../replay.sh) gives the raw extrema and 3,072 checked boxes, zero failures.

## 2. Residual budget

For $v=y,m\ge24$ the [outer proof](outer-moderate.md) gives, on $0\le p\le1/2$,
$$
 |T_v|,|T_v'|\le3(v+2)e^{-15v/16},\qquad |T_v''|\le3(v+2)^2e^{-v},
$$
and for $u_v=\sqrt{W_v}$, the bounds $u_v\ge9/10$, $|u_v'|\le1$,
$|((1+u_v)^{-1})'|\le1$, $|((1+u_v)^{-1})''|\le4$.
The true $u_t\le1$ at all depths follows from $0\le J_t\le e^{-pt}/p$. The [zero-cross argument](zero-cross-helmholtz.md) includes $t=0$ in this statement and in the high-band convexity below.

Use the exact identity
$$
 D=e^{-cp}G+E,\qquad E=-T_m+ Z e^{cp}T_y,\qquad Z=u_x/(1+u_y) .
$$
By (2), $|Z|\le1$, $|Z'|\le4$, $|Z''|\le32$. Put $w=(y+7)/2\ge m+2$, $H=3(y+2)e^{-11y/16}$,
$Q=e^{cp}T_y$. Since $e^{cp}\le e^{y/4}$,
$|Q|\le H$, $|Q'|\le(c+1)H$, $|Q''|\le(c^2+2c+y+2)H$.
Hence $|(ZQ)'|\le(c+5)H$ and $|(ZQ)''|\le(c^2+10c+y+42)H$.
The low residual cost
$ |E'(0)|+\int_0^{1/2}|c^2 E-E''|$
is at most
$$
 3w(1+y^2/8)e^{-15y/32}
 +\tfrac32 w^2 e^{-y/2}
 +3(y+2)(y^2/4+7y/2+26)e^{-11y/16}.
$$
Every summand decreases on $y\ge48$ since its polynomial has nonnegative coefficients, degree at most three, and exponential rate at least $15/32$.
At $y=48$ the sum is $< (45/32)/1000$, hence less than $cg_0/1000$.

Above $\beta=1/2$ set $U=J_m$, $V=\sqrt{J_xJ_y}$. The [true high profiles](outer-moderate.md) are nonnegative, nonincreasing and convex, bounded by $e^{-pm}/p$. At $p=\beta$ one has
$J_t=e^{-pt}u_t^2/p$. Thus $-U'\le(m+5)U$ by $u_m\ge9/10$, $|u_m'|\le1$, and $-V'\le(m+12)V$ by $|u_x'/u_x|\le8$, $|u_y'/u_y|\le10/9$. For a sharp enough common bound on the sum, note that the *first* coefficient is actually $m+2+20/9$, while the second is at most $m+2+8+10/9$; their sum is $2m+46/3<2m+16$. Therefore
$$
 -U'(\beta)-V'(\beta)\le (4m+32)e^{-m/2}.
$$
Monotonicity and convexity then bound $\int_\beta^L|D''|+|D'(L-)|$ by that sum, including the endpoint atom. Also $\int_\beta^L |D|\le2e^{-m/2}/m$. Using $c\ge45/2$, $m/c\le17/15$, $m\ge24$ and (2), the sum of both costs (the integral of $|D|$ multiplied by $c^2$), divided by $cg_0$, is at most
$$
 16[2+4(17/15)+32/(45/2)]e^{-12}<1/1000.
$$

## 3. Conclusion

Extend $G$ affinely past $1/2$. As in [the signed Fourier argument](outer-moderate.md#23-signed-transform), two integrations by parts on the residual and the affine pieces give a baseline $cg_0$, with costs $|G'(0)|$, $\int_0^{1/2}e^{-cp}|G''|\le4/c$, the low and high residuals just bounded, and $2|G'(1/2)|e^{-c/2}$. This last cost is $<cg_0/1000$. Hence the normalized lower bound is at least
$$
 1-\frac{3/5}{45/32}-\frac{64}{(45/2)^2}-\frac3{1000}
 =\frac{35957}{81000}>\frac25.
$$
The same checker independently tests every last numerical fee with rational partial sums for $e^t$ in the denominator. From the project root run:

    PYTHONDONTWRITEBYTECODE=1 python3 certificates/verify_zero_far48.py

The original Fourier sampling then yields the desired finite Hermite comparison on this entire branch, including its zero-depth endpoint.
