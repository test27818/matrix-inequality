# 半正定 Hermitian 矩阵的乘积不等式

**定理。** 对任何正整数 $n$ 和复 Hermitian 半正定矩阵 $A=(a_{ij})_{1\le i,j\le n}$，
$$
\prod_{i,j=1}^n |1-a_{ij}| \le \prod_{i=1}^n |1-a_{ii}^{\,n}|.
$$

这是计算机辅助数学证明。矩阵归约和连续转接是解析的，标量部分使用解析界与可重放的外舍入区间证书；证书覆盖连续参数盒以及完整频轴，不是浮点试点。本文给出完整数学证明；区间程序及其正确性说明与正文一同提供。

## 1. 从标量下界归约到矩阵

$n=1$ 或存在对角元为 1 都直接成立。其余情况令 $d_i=a_{ii}$、$s_i=\max(1,d_i)$，并置 $B_{ij}=a_{ij}/(s_i s_j)$、$r_i=B_{ii}<1$。正对角合同保持半正定。取 $t=s_i s_j,z=a_{ij}$，由 $|z|^2\le d_i d_j\le t$ 得
$$t|1-z/t|^2-|1-z|^2=(t-1)(1-|z|^2/t)\ge0.$$
将各项开方相乘，归约因子是 $\prod_i s_i^n$；右端用 $|1-d_i^n|=s_i^n(1-r_i^n)$ 消去同一因子。再消去两边共有的正因子 $\prod_i(1-r_i)$，目标成为
$$\prod_{i\ne j}|1-B_{ij}|\le\prod_i S_n(r_i),\qquad S_n(r)=\sum_{k=0}^{n-1}r^k.$$

取 $0<r,s<1,x=-n\log r,y=-n\log s$。设 $J_x(p)=\int_x^\infty[\sinh((1-p)v)-(1-p)v]/(\cosh v-1)\,dv$（$0<p<1$）。置
$$g_0(r)=-\frac{x}{n(e^x-1)},\qquad g_j(r)=\frac1n J_x(j/n)\quad(1\le j<n).$$
下面证明对 $|z|\le1$，
$$
-\log|1-\sqrt{rs}z|\ge\frac{g_0(r)+g_0(s)}2+\sum_{j=1}^{n-1}\sqrt{g_j(r)g_j(s)}\Re z^j. \tag{1}
$$
半径为 0 时定义所有 $g_j(0)=0$；(1) 显然仍成立，因为 $g_0\le0$。另有 [常数项恒等式](continuum/constants.md)：
$$(n-1)g_0(r)-\sum_{j=1}^{n-1}g_j(r)=-\log S_n(r). \tag{2}$$
对正半径的行作归一化，写 $B_{ij}=\sqrt{r_ir_j}\,C_{ij}$，其中 $C\succeq0$ 且 $C_{ii}=1$。半径为零的行原来就是零行，在 $C$ 中补对应的单位行列即可。对每个 $j\ge1$，令 $v_i=\sqrt{g_j(r_i)}$；由 Schur 乘积定理，$v^*C^{\circ j}v\ge0$，故非对角项实部之和至少为 $-\sum_i g_j(r_i)$。对所有有序 $i\ne j$ 应用 (1)，按 (2) 合并后指数化，即得所求。

## 2. 为什么连续证书能推出 (1)

记 $m=(x+y)/2$、$D(p)=J_m(p)-\sqrt{J_x(p)J_y(p)}$，在 $p=0$ 连续补值，在 $p\ge1$ 接零。原对数核与极化核之差的连续谱为
$$b(p)=e^{-mp}/p-\mathbf 1_{p<1}\sqrt{J_x(p)J_y(p)}.$$
它是连续可积函数，零点取两深度的正则化值之平均，且在 $p\ge1$ 有指数尾。它是 $D$ 加上同半径谱，后者的余弦变换等于 [显式非负积分](continuum/diagonal.md)：
$$
\frac{1-\cos\omega}{\omega^2}
\int_m^\infty \frac{2v^3e^{-v}}{(1-e^{-v})^2(v^2+\omega^2)}\,dv,
$$
零频率取极限。因此只需证明 $D$ 在任意实频率的余弦变换非负。

把 $b(|t|/n)e^{it\theta}$ 以周期 1 求和。记 $\widehat b(\omega)=\int_0^\infty b(p)\cos(\omega p)\,dp$。周期函数的第 $\ell$ 个傅里叶系数为 $2n\widehat b(n(\theta-2\pi\ell))\ge0$。用 [峰函数积分](continuum/positive-periodization.md) 即得周期函数在 0 的值非负；绝不假定未验证的傅里叶级数收敛。将该值除以 $2n$ 展开，恰是 (1) 的两边之差在 $z=e^{i\theta}$ 上的值。要从圆周推到圆盘，令 $P(z)$ 为 (1) 右侧取实部前的多项式（常数项也包含在内）。对解析函数 $(1-\sqrt{rs}z)e^{P(z)}$ 用最大模原理：它在单位圆周上的模不超过 1，故在整个圆盘内也不超过 1。端点、单位和换极限的细节见 [连续转接](continuum/sampling.md)。

区间证书使用两种充分入口：直接验证上述 $D$（连续的 $\kappa=0$），或验证原有限核 $D_\kappa$，其中密度的分子第二项为 $(1-p)\sinh(\kappa v)/\kappa$，截止为 $1-\kappa$。有限区域证书涵盖 $\kappa$ 趋零；按正则化后的主控收敛，得同一个连续 $D$ 的判号。**不**把连续证书反向当作非零 $\kappa$ 的证书。

## 3. 完整径向覆盖和重放

只考虑 $0<x\le y$；同径 $x=y$ 时 $D=0$。以下每行覆盖全频，交界端点包含在闭盒中。

|较浅深度 $x$|怎样覆盖所有 $y\ge x$|
|---|---|
|$0<x\le3$|[至 24](regions/short-zero-residual.md)，[24 到 48](regions/zero-cross-helmholtz.md)，[48 及以上](regions/zero-far48.md)。|
|$3\le x\le8$|[至 16](regions/compact-continuum-mixed.md)，[16 到 24 的 0、1 号根：从源码重放](../replay.sh)，[24 到 32](regions/compact-cross32.md)，[32 到 48](regions/shallow-cross-helmholtz.md)，[48 及以上](regions/far-shallow-bridge-proof.md)。|
|$8\le x\le12$|[至 16](regions/compact-continuum-mixed.md)，[16 到 24](regions/compact-balanced-proof.md)，上述 24 到 48，加 [48 及以上](regions/far-shallow-physical-proof.md)。|
|$12\le x\le16$|[两深度之差至多 12](regions/compact-near-spectrum.md)，已包含 $y\le24$；其余同上一行。|
|$16\le x\le24$|[间距至多 8](regions/compact-near-spectrum.md)，[间距 8 到 16](regions/compact-gap-helmholtz.md)，两者已包含 $y\le32$；加 [32 到 56](regions/moderate-helmholtz.md) 和 [56 及以上](regions/outer-moderate.md)。|
|$x\ge24$|[间距至多 18](regions/deep-continuum-closure.md)（压缩无界深度），[间距至少 18](regions/deep-separation-proof.md)。|

用 [精确覆盖脚本](../certificates/verify_cover_geometry.py) 可独立查验这些区域的并确实是整个正半径三角域。证明性程序入口：

```sh
sh replay.sh all
```

它重建各需要的 C++20 可执行文件，并调用 Python 3 标准库的精确有理/整数检查。环境需 IEEE binary64（禁快速浮点和乘加融合）、GMP 头文件和 MPFR 6 共享库；具体外舍入与误差界见每个分区的说明。旧日志、未完成的探索和带额外前提的 Lean 文件都不是这里的证明前提。
