# 远浅半带的物理标量不等式：不依赖旧区间树的完整证明

**定理。** 设 $8\le x\le16,\ y\ge48,\ 0\le\kappa\le1/16$，并记 $m=(x+y)/2,\ c=(y-x)/2,\ a=y/2,\ L=1-\kappa,\ \beta=1/2$。令 $S_k(v)=\sinh(kv)/k$（$S_0(v)=v$），并在 $0<p<L$ 定义
$$
 J_v(p)=\int_v^\infty\frac{\sinh((1-p)t)-(1-p) S_\kappa(t)}{\cosh t-1}\,dt,
 \qquad D(p)=J_m(p)-\sqrt{J_x(p)J_y(p)} .
$$
在 $p=0$ 对 $D$ 取右极限，$p\ge L$ 接零。置
$g_0=\tfrac12\lim_{p\downarrow0}(e^{-px}/p-J_x(p))$。则 $g_0\ge xe^{-x}>0$，且对**每个实频率** $\omega$，
$$
 (c^2+\omega^2)\int_0^\infty D(p)\cos(\omega p)\,dp
 \ \ge\ \frac{24167}{160000}\,c g_0\ >\ \frac18 c g_0. \tag{1}
$$
这里 $y$ 没有上限，频率和 $\kappa$ 不需离散采样。这是高阶矩阵证明用到的一个区域，不是整个矩阵定理。以下全是实分析与有理常数估计；最后的标准库 Python 检查器只核算常数，不充当分析引理。

## 1. 正则级数、去极点和精确分解

写 $\mu=1-p$。由于 $S_\mu(t)$ 对 $\mu\ge0$ 非减，分子 $\mu(S_\mu-S_\kappa)$ 非负。密度还不超过 $e^{-pt}$：将差乘以 $\cosh t-1$，得到
$$
 \mu S_\kappa(t)+e^{-t}\cosh(pt)-e^{-pt}\ge0 .
$$
确实 $e^{-pt}-e^{-t}=e^{-pt}(1-e^{-\mu t})\le\mu t$。所以 $0\le J_v(p)\le e^{-pv}/p$。

设 $q=e^{-v}, h=\kappa^2, S=S_\kappa(v), C=\cosh(\kappa v)$。正则化 $T_v=e^{-pv}/p-J_v$ 在 $p=0$ 亦有下列展开（下标 $v$ 暂省）：
$$
 T=(1-p)A+B-E_*,\quad A=\sum_{j\ge1}\frac{2j q^j(jS+C)}{j^2-h},
 \quad B=\sum_{j\ge1}\frac{j q^j e^{-(1-p)v}}{j+1-p},
 \quad E_*=\sum_{j\ge1}\frac{(j+1)q^j e^{-pv}}{j+p}. \tag{2}
$$
证明：展开 $(\cosh t-1)^{-1}=2\sum_{j\ge1}j e^{-jt}$；其中 $e^{\mu t}$ 列的 $j=1$ 项恰为 $e^{-pv}/p$，其余同列换编号 $j\mapsto j+1$。另有
$\int_v^\infty e^{-jt} S_\kappa(t)dt=q^j(jS+C)/(j^2-h)$，在 $\kappa=0$ 也成立。所有积分在 $p>0$ 绝对收敛。在 $[0,\beta]$ 稍扩大的开区间上，剩余级数及其一二阶导数均受固定倍数几何级数控制（$v$ 固定）。$B,E_*$ 的一阶导数分别逐项乘以 $v+1/d,-(v+1/d)$，二阶导数均乘 $(v+1/d)^2+d^{-2}$，$d$ 表示相应分母 $j+1-p$ 或 $j+p$。

定义 $W_v=1-p e^{pv} T_v$；$p>0$ 时 $W_v=p e^{pv}J_v\in[0,1]$，$W_v(0)=1$。在低段 $[0,\beta]$ 有
$$
 D=e^{-cp}G+E,\quad
 G=\frac{T_x}{1+\sqrt{W_x}},\quad
 E=-T_m+K e^{cp}T_y,\quad
 K=\frac{\sqrt{W_x}}{1+\sqrt{W_y}} . \tag{3}
$$
事实上 $D=e^{-pm}(1-\sqrt{W_xW_y})/p-T_m$，用
$1-\sqrt{W_xW_y}=(1-\sqrt{W_x})+\sqrt{W_x}(1-\sqrt{W_y})$ 得 (3)，包括 $p\downarrow0$ 的极限。于是 $G(0)=T_x(0)/2=g_0$。

## 2. 浅端主项：零模首坡和全低段曲率

本节取 $v=x, u=v e^{-v}$。由 $0\le\kappa v\le1$、$\sinh z/z=\int_0^1\cosh(sz)ds$，
$v\le S\le1.18v,\ 1\le C\le1.55,\ q<.00036,\ e^{-v/2}<.019$。逐列加几何级数：
$$
 A\le\frac{2q(S+C)}{(1-h)(1-q)}\le2.77u,\quad
 B\le .02q,\quad E_*\le2.01q .
$$
$T'=-A+B'+|E_*'|$，$T''=B''-E_*''$。因而
$$
 |T|,|T'|\le2.81u,\quad |T''|\le2.01q(v^2+2v+2)\le37u .
$$
最后用了 $(v^2+2v+2)/v\le18.25$。在 $p=0$，每个 $A_j-E_{*,j}\ge q^j(2S+1/j-1)\ge0$，第 1 项至少 $2Sq$，故 $2Sq\le T(0)\le2.81u\le1$、$g_0\ge u$。

如果 $\kappa=0$，令 $A^0=A(0,v,1)$；在 $p=0$，
$$
 A^0_j-|E_{*,j}'|=q^j(1-1/j)(v+1/j).
$$
$j=1$ 时为零，$j\ge2$ 时恰为 $B'_{j-1}$。所以任意 $\kappa$ 都有 $T'(0)=-(A-A^0)$。每个差除以 $q^j$ 后为
$$
 \frac{2}{1-h/j^2}\left[S-v+(C-1)/j+(h/j^2)(v+1/j)\right]\in[0,S(1-q)/2].
$$
上界用 $S-v\le9S/59, C-1\le11S/160, h(v+1)\le9S/(8\cdot256)$ 代入即可。求和得 $-Sq/2\le T'(0)\le0$。由 $G'(0)=T'(0)/2+T(0)^2/8$ 得
$$
 |G'(0)|\le g_0/4. \tag{4}
$$

设 $Z=p e^{pv}T$，则
$$
 Z'=e^{pv}[(1+pv)T+pT'],\qquad
 Z''=e^{pv}[(2v+pv^2)T+(2+2pv)T'+pT''] .
$$
用 $v^j e^{-v/2}\le8^j e^{-4}$（$1\le j\le3,v\ge8$），上述三式给
$(|Z|,|Z'|,|Z''|)\le(.22,2.4,28)$，例如最后一个界来自
$[2.81(2v+v^2/2+2+v)+37/2]v e^{-v/2}$。现在 $r=\sqrt W\ge.85$。
$\ell(w)=(1+\sqrt w)^{-1}$ 的函数值及前两阶导数的绝对值上界为 $.55,.18,.25$：其中 $|\ell'|=[2r(1+r)^2]^{-1}$，
$\ell''=[4r^3(1+r)^2]^{-1}+[2r^2(1+r)^3]^{-1}$。代入 $G=T\ell(W)$ 得
$$
 |G''|/u\le37(.55)+2(2.81)(.18)(2.4)+2.81[(.18)28+(.25)2.4^2]<42 .
$$
特别地
$$
 |G''|\le45g_0,\qquad |G'(\beta)|\le(1/4+45/2)g_0. \tag{5}
$$

## 3. 高段的变差和截止原子

为避免在 $p=L$ 时对零除法，先除去显式因子 $P(\mu)=\mu(\mu^2-\kappa^2)$。有
$$
 \sinh(\mu t)-\mu S_\kappa(t)
 =P(\mu)\int_0^t S_\mu(t-s)S_\kappa(s)ds. \tag{6}
$$
以 $f(s)=S_\mu(t-s),g(s)=S_\kappa(s)$ 积分 $(f'g-fg')'=(\mu^2-\kappa^2)fg$ 即得（当两参数相等也成立）。将卷积除以 $\cosh t-1$ 再从 $v$ 积至无穷，记为 $A_v^+(\mu)>0$。于是 $J_v=P A_v^+$。在 $\kappa\le\mu\le1/2$ 的邻域它光滑：在 $|\mu|,\kappa\le r<1$，对 $\mu$ 求 $j$ 次导的被积函数受 $C_{v,j}(1+t)^{j+3}e^{-(1-r)t}$ 控制。由 $S_\mu(s)=\frac12\int_{-s}^s e^{\mu w}dw$，$A_v^+$ 为与 $\mu$ 无关的正测度的指数矩，故由 Hölder 不等式对数凸；也因 $S_\mu(s)=\int_0^s\cosh(\mu w)dw$ 而对 $\mu\ge0$ 非减。其几何均值 $\sqrt{A_x^+A_y^+}$ 同样非减且对数凸，因而凸。$P=\mu(\mu-\kappa)(\mu+\kappa)$ 在本段非负、非减、凸。非负非减凸函数之积仍凸。故
$$
 U=P A_m^+=J_m,\qquad V=P\sqrt{A_x^+A_y^+}=\sqrt{J_xJ_y}
$$
关于 $p=1-\mu$ 均凸而非增，直到 $L$ 的单侧导数都是有限的。这里**不需要** $D=U-V$ 凸。于是
$$
 \int_\beta^L |D''|dp+|D'(L-)|\le
 \int_\beta^L(U''+V'')dp-U'(L-)-V'(L-)=-U'(\beta)-V'(\beta). \tag{7}
$$

下面给出真半点首坡，无需调用旧数值树。对每个 $v\ge8$ 写 $z=e^{v/16},\ d=e^{-v/2}$。$S\le v z, C\le z$，故 (2) 给
$A\le2.04(v+1)zq$。在 $\beta$，$B\le1.01 d q,E_*\le2.01d q$。保留 $-A$ 的符号，得
$$
 T(\beta)\le1.2(v+1)zq,\qquad T'(\beta)\le3.1(v+1)d q.
$$
从而 $T'(\beta)+(v+6)T(\beta)\le8d$：除以 $d$ 后上界是
$3.1(v+1)e^{-v}+1.2(v+6)(v+1)e^{-7v/16}$。利用
$v e^{-v}\le8(.019)$、$v^2e^{-7v/16}\le64 e^{-7/2}<64/23$、
$v+1\le9v/8, (v+6)(v+1)\le126v^2/64$，上界不超过 $3.1\cdot9\cdot.019+1.2\cdot126/23<8$。因为 $J_v=e^{-pv}/p-T_v$，这正是
$-J_v'(\beta)\le(v+6)J_v(\beta)$。给几何均值求对数导（在 $\beta$ 各 $J_v>0$）后，用 $U(\beta),V(\beta)\le2e^{-m/2}$ 得
$$
 -U'(\beta)-V'(\beta)\le(4m+24)e^{-m/2}.
$$
此外 $0\le U,V\le e^{-pm}/p$，故高段 $|D|\le2e^{-pm}$。记
$$
 H=\frac{c^2\int_\beta^L|D|dp+\int_\beta^L|D''|dp+|D'(L-)|}{c g_0}.
$$
由以上估计和 $g_0\ge x e^{-x}$，
$$
 H\le\left(2c/m+4m/c+24/c\right)e^{(x-c)/2}/x\le21/32 . \tag{8}
$$
最后一步中 $c\ge16,1\le t=m/c=1+x/c\le2$，所以 $2/t+4t\le9$。又 $c\ge24-x/2$，且
$16/x=1+(16-x)/x\le e^{3(16-x)/4}$，故指数商不超过 $1/16$。

## 4. 全低段的深端余项

对所有 $v\ge8,\ 0\le p\le\beta$，由 (2)，
$A\le2.04(v+1)e^{v/16}q,\ B\le .02q,\ E_*\le2.01q$。按第 2 节的一阶符号比较及二阶的正列相减，可得
$$
 |T_v|,|T'_v|\le3(v+2)e^{-15v/16},\qquad
 |T''_v|\le3(v+2)^2e^{-v}. \tag{9}
$$
例如 $B'+|E_*'|\le2.03(v+1)q$，二阶则由 $2.01(v^2+2v+2)q$ 控制。对 $v\ge32$ 写 $s=e^{-7v/16}$，乘 $e^{pv}$ 后 (9) 三界可取 $6vs,6vs,12v^2s$。代入第 2 节的 $Z,Z',Z''$ 公式，得
$$
 |Z|\le3vs,\quad |Z'|\le6v(1.5+v/2)s,\quad
 |Z''|\le[6v(2v+v^2/2)+6v(2+v)+6v^2]s .
$$
对 $1\le j\le3$，$v^j e^{-7v/16}\le32^j e^{-14}<32^j/10^6$，故三式分别不超过 $.22,1,1$。结合第 2 节的浅端估计，$R=\sqrt{W_x}$ 的三阶位（函数、一阶、二阶）绝对值不超过 $(1,2,22)$，其中 $|R''|\le28/(2\cdot.85)+2.4^2/(4\cdot.85^3)$。深端 $\ell(W_y)$ 的对应三界可取 $(1,1,1)$。因而 $K=R\ell(W_y)$ 的三界不超过 $(1,3,28)$。

设 $Y=y+2$，$Q=e^{cp}T_y,\ g=e^{-11y/16-x/4}$。$c\le Y/2$ 且 $e^{cp}e^{-15y/16}\le g$，故
$$
 (|Q|,|Q'|,|Q''|)\le(3Y,2Y^2,2Y^3)g .
$$
例如二阶由 $(3c^2Y+6cY+3Y^2)g$ 控制，且 $Y\ge50$。写 $A_0=a+10$；因 $m+2\le A_0$，将 (9) 用于 $v=m$。再用 $g_0\ge x e^{-x}$、$e^{17x/32}/x<700,\ e^{x/2}/x<400,\ e^{3x/4}/x<22000$（用 $x\ge8,x\le16$ 分别界定分母、分子即足够）。对 $E=-T_m+KQ$ 的三阶位有 $|E^{(j)}|/g_0\le r_j$，其中
$$
\begin{aligned}
 r_0&=3A_0\,700 e^{-15a/16}+3Y\,22000 e^{-11a/8},\\
 r_1&=3A_0\,700 e^{-15a/16}+(9Y+2Y^2)22000 e^{-11a/8},\\
 r_2&=3A_0^2\,400 e^{-a}+(84Y+12Y^2+2Y^3)22000 e^{-11a/8}.
\end{aligned}
$$
由于 $2a/3\le c\le a$，
$$
 B_E:=\frac{|E'(0)|+\int_0^\beta|c^2E-E''|dp}{c g_0}
 \le\frac a2r_0+\frac{3r_1}{2a}+\frac{3r_2}{4a}<\frac3{10000}. \tag{10}
$$
严格上界并非抽样。置 $t=a/24\ge1$，则 $A_0\le34t,Y\le50t$；三个指数前的多项式及有理因子分别至多为 $t^3$ 乘在 $a=24$ 的系数（分母里的 $t$ 可丢弃）。对于 $b\ge1/8$ 有 $t^3 e^{-ba}\le e^{-24b}$。在三个指数处把该上界依次放大到 $(10/27)^{22},(10/27)^{24},(10/27)^{33}$，利用 $e>27/10$。所得有理数为
$$
 \frac{48451154913004447727875000000000000000000000}
 {171792506910670443678820376588540424234035840667}<3/10000 .
$$

## 5. Fourier 合并：所有频率同时成立

$G_*$ 在 $[0,\beta]$ 取 $G$，其后按切线 $G(\beta)+(p-\beta)G'(\beta)$ 延拓。于是
$$
 G_*(p)=g_0+pG'(0)+\int_0^\beta(p-s)_+G''(s)ds .
$$
记 $I(f)=\int_0^\infty f(p)\cos(\omega p)dp$。用 $z=c-i\omega$ 以及 $\int_0^\infty e^{-zp}dp=1/z,\ \int_0^\infty p e^{-zp}dp=1/z^2$，交换绝对可积项后取实部，可得
$$
 (c^2+\omega^2) I(e^{-cp}G_*)\ge c g_0-|G'(0)|-\int_0^\beta e^{-cp}|G''(p)|dp .
$$
再令 $R=D-e^{-cp}G_*$，在低段它等于 $E$；在 $\beta$ 一阶相容（实际 $D$ 在此光滑）。$D$ 在 $L$ 连续、右侧恒零，导数跳跃是 $-D'(L-)$。因此两次逐段分部积分得
$$
 (c^2+\omega^2) I(R)=-R'(0)+\int_0^\infty(c^2R-R''_{\mathrm{pw}})\cos(\omega p)dp+D'(L-)\cos(\omega L).
$$
$R''_{\mathrm{pw}}$ 表示开区间上通常的二阶导，未暗中丢失截止原子。高段及截止以后的切线延拓项满足
$(c^2-\partial_p^2)(e^{-cp}G_*)=2c e^{-cp}G'(\beta)$，绝对值的积分是 $2e^{-c/2}|G'(\beta)|$。用 (4)、(5)、(8)、(10) 和 $c\ge16$，就得到
$$
 \frac{(c^2+\omega^2) I(D)}{c g_0}\ge
 1-\frac1{64}-\frac{45}{256}-\frac1{1000}-\frac3{10000}-\frac{21}{32}
 =\frac{24167}{160000}>\frac18 .
$$
其中切线费 $2e^{-c/2}|G'(\beta)|/(cg_0)\le2e^{-8}(1/4+45/2)/16<1/1000$。这证明了 (1)，包括 $\kappa=0$。没有使用黑箱区间证书或假设真实谱的导数估计。有限阶下再取 $\kappa=1/N$ 做 Fourier 取样；取样、Hermite 系数识别与矩阵张量平方求和是 (1) **之后**的步骤，不是 (1) 的前提。

## 复核

执行 `python3 certificates/verify_far_shallow_physical.py`。只使用 Python 标准库的分数。对 $x\ge0$ 的指数以 Taylor 下和与首个余项乘以 $1/(1-x/(n+2))$ 夹逼，$x<0$ 则交换端点取倒数。输出应含 `normalized guaranteed margin = 24167/160000 > 1/8`。已分块通过 Lean 的主要引理见根目录 ShallowDerivative、ShallowSlope、ZeroMode、SharpOrdinate、PositiveConvolution、HighSpectrum、HighProfiles、HalfSlope、HalfProfiles、DeepTail、DeepJet 和 ConeFactors。本文是独立成篇的**数学证明**；不声称有一条无条件 Lean 定理已一键串起整个 (1) 或全阶矩阵命题。
