# 中段闭窗：保留高谱差与端点的独立精确核验

在 $3\le x\le16$, $16\le y\le24$, $x<y$, $0\le b=\kappa^2\le1/256$，令
$c=(y-x)/2$、$m=(y+x)/2$、$L=1-\kappa$。本页验证**有限阶真实谱**
$D=J_m-\sqrt{J_xJ_y}$ 在 $L$ 后置零的余弦变换，连接见[采样文档](../continuum/positive-periodization.md)。只有下表实际通过的整体闭窗计入覆盖，不计被拒绝的其他根或其中个别叶子。

## 去极点、截止与全相位判号

小 $p$ 使用[中深文档公式 (1)](moderate-helmholtz.md)：
$$
D=e^{-cp}G+E,\quad G=T_x/(1+\sqrt{W_x}),\quad
E=-T_m+e^{cp}\sqrt{W_x}T_y/(1+\sqrt{W_y}).
$$
$T_t=e^{-pt}/p-J_t$ 用去极点级数算。一般乘子 $\lambda>0$ 下以
$$
D''-\lambda^2D=e^{-cp}\{(c^2-\lambda^2)G-2cG'+G''\}+E''-\lambda^2E
$$
作区间计算，$\lambda=c$ 时第一个差额严格为零，不必对两个区间相减。

高 $p$ 不分别支付 $J_m$ 与 $\sqrt{J_xJ_y}$ 的总变化。置 $q=1-p$，$P=q(q^2-b)$，
$$D=P(A_m-\sqrt{A_xA_y}),\qquad
A_t=2\sum_{n\ge0}\sum_{l\ge1}l M_{2n+3}(t,l)\sum_{i=0}^n q^{2i}b^{n-i},\quad
M_j(t,l)=\int_t^\infty \frac{v^j}{j!}e^{-lv}\,dv.$$
这是[原谱正卷积](moderate-helmholtz.md)的正泰勒展开，根号分母在 $q=\kappa$ 不丢失正性。
延拓区 $q<\kappa$ 函数及其普通二阶导均为零；区间若可能跨越截止点，取包含零的下界。
额外的分布原子 $-D'(L-)$ 始终非负，这是[独立的全深度解析引理](cutoff-atom-logconcavity.md)；
Python 简单预算还显式扣除数值区间可能给出的不利原子界，因而它本身不依赖这个引理。

分部积分得，对于 $\omega\in\mathbb R$，
$$
(\lambda^2+\omega^2)\int_0^L D(p)\cos(\omega p)\,dp
=\lambda^2\int_0^L D(p)\,dp+\int_0^L(D''-\lambda^2D)(1-\cos(\omega p))\,dp
-D'(L-)(1-\cos(\omega L)).
$$
在连续谱闭格上包住 $\lambda^2D$、$D''-\lambda^2D$ 的下端，给出以下全频证书。
* 若对所有角度共用 $0\le1-\cos\le2$ 已足够，立即取积分下和。Python 版本直接取同一点两个区间下端 $(\lambda^2D,2D''-\lambda^2D)$ 的较小者，以保留变量依赖。
* 否则令 $d_i,h_i$ 分别为两项在第 $i$ 格的下界，以 $M=\sum(d_i+h_i)\Delta_i$、$V=|h_{K-1}|+\sum_{i=1}^{K-1}|h_i-h_{i-1}|$。高频由逐格积分望远镜求和，下界至少 $M-V/|\omega|$。其余区间用宽度 $1/2$ 的闭频率格，积分权重在频率中心 $w$、谱格中心 $z$、半宽 $r$ 时恰为 $1-\cos(wz)\operatorname{sinc}(wr)$；整个频率格误差至多 $z/4$。每个参数闭叶可各自选择 $\lambda$ 或某个证法，原 $D$ 则始终相同。

## 两个算术实现

**C++ 二进区间相位核。** [compact_angle_residual.cpp](../../certificates/compact_angle_residual.cpp) 直接复用 [short_zero_residual.cpp](../../certificates/short_zero_residual.cpp) 的正矩、极点消除和算术。低频把 $l=1$ 整项精确抽出，$l\ge2$ 的余项在 $q\le1$ 也按 $1/4$ 的几何率衰减；$q\le1/2$ 时全部 alias 的 $q$ 级数有同样几何尾。对重叠频段同时计算两种外包，取交缩紧。径向端点、频谱格和切点可用二进制精确表示；四则用 IEC559 双精度，每次向外扩一相邻表示数，禁 fast-math/FMA；平方根候选再以两次外包平方复核，指数和三角函数均用带余项的 Taylor 和二倍角恢复。编译启用 UBSan，无通用 libm 的三角/指数近似作为证明。依赖请保留同目录两个源文件。

## 用于主定理的闭窗

$[8,12]	imes[16,24]$，全部 $\kappa^2\in[0,1/256]$ 及实频率。
651 次访问、326 个接受闭叶，其中 12 叶用直接符号，314 叶用连续角度区间。

```sh
cd certificates
c++ -O3 -std=c++20 -ffp-contract=off -fno-fast-math -fsanitize=undefined -fno-sanitize-recover=all compact_angle_residual.cpp -o /tmp/compact-angle
/tmp/compact-angle 8 12 16 24 32 24 80000
```
