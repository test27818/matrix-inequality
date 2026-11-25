# 中深交叉窗的独立 Helmholtz 重建

**结论。** 对任意 $b=\kappa^2\in[0,1/256]$，下列三个闭盒使用同一真实径向谱，均以严格正余量证明中点 Hermite 多项式支配几何均值极化多项式：

* $16\le x\le20$, $32\le y\le48$；
* $20\le x\le24$, $32\le y\le48$；
* $16\le x\le24$, $48\le y\le56$。

这是新写的解析降维加可重放整数区间证明，不是找回的旧 C++ 树。这里的 Python 重放不能代替 Lean 内的真实谱定理；独立的 Lean 重放版本及其语义连接，以工程总说明为准。只覆盖这三个中深交叉窗，不自动推出其他区域。

## 1. 真谱、消极点、截止

置 $c=(y-x)/2$, $m=(x+y)/2$, $L=1-\kappa$, $\mu=1-p$，并令 $G_a(v)=\sinh(av)/a$，在 $a=0$ 取 $G_0(v)=v$。在 $0<p<L$ 定义
$$
 F_p(v)=\frac{\mu(G_\mu(v)-G_\kappa(v))}{\cosh v-1},\quad J_t(p)=\int_t^\infty F_p(v)\,dv,\quad
 D(p)=J_m(p)-\sqrt{J_x(p)J_y(p)}.
$$
取 $E_t(a)=e^{-at}/a$，只在 $p>0$ 时先写 $T_t(p)=E_t(p)-J_t(p)$。令 $W_t=1-p e^{pt}T_t$。由下节的绝对收敛级数，$T_t$ 可延至 $p=0$，且 $W_t(0)=1$。
对 $0<p<L$ 的精确恒等式是
$$
 D=e^{-cp}\frac{T_x}{1+\sqrt{W_x}}-T_m+
 e^{cp}\frac{\sqrt{W_x}}{1+\sqrt{W_y}}T_y.                         \tag{1}
$$
因此 $D(0)=(T_x(0)+T_y(0))/2-T_m(0)$；同一式子给零频率处的连续一、二导数，不对原来含 $1/p$ 的式子作区间运算。

高频端记 $P=\mu(\mu^2-b)$。利用 $G_a''=a^2G_a$, $G_a(0)=0$, $G_a'(0)=1$，卷积的初值和二阶导数给
$$
 \frac{G_\mu(v)-G_\kappa(v)}{\mu^2-\kappa^2}
   =\int_0^v G_\mu(v-s)G_\kappa(s)\,ds .                         \tag{2}
$$
右式还给 $\mu=\kappa$ 的连续极限。故 $A_t=J_t/P$ 在 $p=L$ 有正的光滑延拓；$D=P(A_m-\sqrt{A_xA_y})$。特别地 $D(L)=0$，$D'(L-)=-2b(A_m-\sqrt{A_xA_y})\bigr|_L$，在 $b=0$ 该斜率为零。以下让 $D(p)=0$ 对 $p>L$，不会误把截止的一阶导数跳跃删掉。

可供独立核对 $A_t$ 的显式正系数式为
$$
 A_t(p)= 2\sum_{q\ge1}q\sum_{j\ge1} M_{2j+1}(t,q)
        \sum_{h=0}^{j-1}\mu^{2h}b^{j-1-h},\quad
 M_\ell(t,a)=\int_t^\infty \frac{v^\ell}{\ell!}e^{-av}dv
 =\frac{e^{-at}}{a^{\ell+1}}\sum_{i=0}^{\ell}\frac{(at)^i}{i!}.
$$
这里正项允许逐项积分；该式不是当前脚本需要使用的高频计算。式 (2) 更直接给下一节需要的凸性。

## 2. 无角度的余量和高频端的统一预算

令 $H=c^2D-D''$，再把 $D'(L-)$ 作为 $p=L$ 的原子。两次分部积分，对每个 $\omega\in\mathbb R$ 得
$$
 (c^2+\omega^2)\int_0^L D(p)\cos(\omega p)\,dp
 =-D'(0)+\int_0^L H(p)\cos(\omega p)\,dp+D'(L-)\cos(\omega L).
$$
其下界即**符号测度余量**
$$
 \int_0^L \big[c^2D-2\max(H,0)\big]\,dp-2\max(D'(L-),0).          \tag{3}
$$
这是因为原子加密度的总质量等于 $c^2\int D+D'(0)$，任一余弦值介于 $-1,1$。不用假设 $D''$、$H$ 或截止原子的号。

取固定分点 $\beta=29/32<L$。在 $\beta\le p\le L$，式 (2) 将 $A_t$ 写成 $\mu$ 的正指数矩：$G_\mu(s)=\frac12\int_{-s}^s e^{\mu u}\,du$。径向积分及一、二导数受 $O(v^5 e^{-(1-3/32)v})$ 控制，可以微分。对称性给 $A_t$ 随 $\mu\ge0$ 非减，Cauchy–Schwarz 给对数凸性；$\sqrt{A_xA_y}$ 同样。$P,P_\mu,P_{\mu\mu}$ 在 $\mu\ge\kappa$ 都非负。因此 $U=J_m$ 和 $V=\sqrt{J_xJ_y}$ 关于 $p$ 均凸且不增；它们在 $L$ 处取零并有左导数。
在此高频段还有 $0\le F_p(v)\le e^{-pv}$：上界改写为
$$
 (2-e^{-v})e^{-pv}\le e^{-\mu v}+2\mu G_\kappa(v).
$$
左端 $\le2e^{-\beta v}\le1$（$v\ge16$, $\beta\cdot16\ge1$, $e^1>2$）；右端 $\ge (1-\mu v)+2\mu v\ge1$。故 $U,V\in[0,e^{-pm}/p]$，$|D|\le e^{-pm}/p$。逐点
$$
 c^2D-2\max(H,0)=\min(c^2D,2D''-c^2D)\ge-c^2|D|-2V''.
$$
又 $\max(D'(L-),0)\le -V'(L-)$。于是 (3) 后半段连同截止原子整体不小于
$$
 2V'(\beta)-\frac{c^2 e^{-\beta m}}{\beta m}.                       \tag{4}
$$
这是对整条高频段和原子的联合下界，不是只检查频率截点。在 (4) 中只须用 (1) 的同一 $W$ 在固定点 $p=\beta$ 计算 $V'$：
$$
 V(p)=\frac{e^{-pm}}p\sqrt{W_x(p)W_y(p)}.
$$

## 3. 低频段的显式有理检查式

恒等式 $(\cosh v-1)^{-1}=2e^{-v}\sum_{q\ge1}q e^{-(q-1)v}$ 给
$$
 T_t(p)=\mu H_t+\sum_{q\ge1}q E_t(q+1-p)-\sum_{q\ge1}(q+1)E_t(q+p),
 \quad H_t=2\sum_{q\ge1}q\int_t^\infty e^{-qv}G_\kappa(v)\,dv.      \tag{5}
$$
所有项都远离 $p=0$ 极点。对 $p\in[0,\beta]$，扣掉 $E_t(p)$ 后，两项移动指数的慢速率至少 1，$G_\kappa$ 项至少 $1-1/16$；在该范围逐项微分两次由绝对一致收敛保证。右侧同时定义 $T_t(0)$。

程序以 $e^t$ 缩放 (5)，二阶 jet 使用 $(f,f',f''/2)$。记 $\widehat E_i(t,a)=e^tM_i(t,a)$。它的一、二号槽在 $a=q+1-p$ 取正号，在 $a=q+p$ 取 $(-1)^i$。两组尾的每个**绝对值**槽分别不超过
$$
 \delta\widehat E_i(t,2-p),\quad \delta\widehat E_i(t,1+p),\qquad
 \delta=\frac{3e^{-t}}{(1-e^{-t})^2}.
$$
理由：每个积分指数差 $q-1$ 对 $v\ge t$ 给乘子 $e^{-(q-1)t}$，而 $\sum_{r\ge1}(r+2)e^{-rt}\le\delta$。分别按原符号把尾放进对应槽。

对于 $e^t H_t$，仅保留 $q=1$ 主项；$q\ge2$ 的正尾至多为完整主项的 $2e^{-t}/(1-e^{-t})^2$ 倍，同理。主项用 $G_\kappa(v)=\sum_{j\ge0}b^jv^{2j+1}/(2j+1)!$ 截至 $j=23$。在本域 $16\le t\le56$，用 $M_i(t,1)\le M_i(0,1)=1$, $e<3$，尾界严格为
$$
 2 e^{56}\frac{(1/256)^{24}}{1-1/256}
 <2\cdot3^{56}\frac{(1/256)^{24}}{1-1/256}<2^{-100}.
$$
右端比较可由整数直接交叉相乘。主项第 $j$ 个积分的指数缩放值为 $\sum_{i=0}^{2j+1}t^i/i!$；不用在 $\kappa=0$ 除以零。

为避免差分丢精度，令 $\widehat T_t=e^t T_t$, $\widehat G=\widehat T_x/(1+\sqrt{W_x})$，$\widehat E=-e^{-c}\widehat T_m+e^{c(p-2)}\sqrt{W_x}\widehat T_y/(1+\sqrt{W_y})$。则
$$
 e^xD=e^{-cp}\widehat G+\widehat E,\qquad
 e^xH=e^{-cp}(2c\widehat G'-\widehat G'')+c^2\widehat E-\widehat E'' . \tag{6}
$$
在每个单元 $[j\beta/K,(j+1)\beta/K]$ 上，程序以经区间延拓的 (5)、(6) 同时包围该单元全体 $p$ 及当前参数盒全体 $(x,y,b)$，不是在端点或中点采样。用 $\min(c^2 e^xD,c^2 e^xD-2e^x H)$ 的区间下端乘长度作为 Riemann 下和项。最后加 (4) 乘 $e^x$ 的区间下端。只有和严格大于零才接受。

## 4. 实现、重放与精确状态

实现：[moderate_helmholtz.py](../../certificates/moderate_helmholtz.py)。所有计算端点为 $2^{-112}$ 的整数倍。加减、乘除直接做整数外舍入；非负平方根用整数平方根；正指数缩小至 $[0,1]$ 后用到 32 次的 Taylor 下界与项比界 $x/34$ 控制余尾，再逐次区间平方，负指数取区间倒数。平方根二阶槽或除法若无法证其正分母，则该参数盒不会通过；程序递归尝试精确二分。每次二分两个孩子共享闭中点，覆盖父盒。所有网格端点和长度精确。

```sh
python3 certificates/moderate_helmholtz.py --root all --mesh 6
```

| 区域 | 参数树访问 | 叶 | 最大深度 | 最小正余量分子，公分母 $2^{112}$ |
|:--|--:|--:|--:|--:|
| $[16,20]\times[32,48]$ | 111 | 56 | 7 | 205018695430871142349943416143930632 |
| $[20,24]\times[32,48]$ | 49 | 25 | 11 | 2562207715038247799653203218317863 |
| $[16,24]\times[48,56]$ | 63 | 32 | 7 | 780363019522303984185032842771470310 |

这些数是 (3) **经过正因子 $e^x$ 缩放**之后的下界。重放日志：[moderate-helmholtz-replay.txt：从源码重放](../../replay.sh)。它们不应与旧代码的权重、分块或余量分子直接比较。严格正余额具体给出 $e^x(c^2+\omega^2)\int D\cos(\omega p)>49/100$，对上述所有参数成立。


## 5. 对接离散 Hermite 核

取偶延拓的 $D$，式 (3) 与分部积分也给 Fourier 变换 $O((1+\omega^2)^{-1})$。由 Fourier 反演，采样序列 $D(j/N)$ 正定。更具体地，对 $M$ 个等比系数作二次型并除以 $M$，得到非负的
$$
 \sum_{|j|<M}(1-|j|/M)D(j/N)e^{ij\theta};
$$
令 $M\to\infty$，因 $D$ 紧支撑，得到所需 $D(0)+2\sum_{j=1}^{N-2}D(j/N)\cos(j\theta)\ge0$。这一步不要求角度采样或原代码的周期化实现。

在 $\kappa=1/N$ 时，同半径 Hermite 多项式的第 $j\ge1$ 模为 $h_j(e^{-t/N})=J_t(j/N)/N$，常数为 $h_0(e^{-t/N})=-T_t(0)/(2N)$。常数可从接触归一化直接校验：由 (5) 或在正 $p$ 微分后取极限得 $\partial_t T_t(0)=F_0(t)-1$；对 $r=e^{-t/N}$，直接求有限几何和，
$$
 \sum_{j=1}^{N-1} F_{j/N}(t)-\tfrac{N-1}{2}\big(F_0(t)-1\big)
     =\frac{r}{1-r}-\frac{N r^N}{1-r^N},\qquad F_{(N-1)/N}(t)=0.
$$
因此 $(N-1)(-T_t(0)/(2N))-\sum_{j=1}^{N-2} J_t(j/N)/N$ 的深度导数与 $-\log(1+r+\cdots+r^{N-1})$ 相同；两边在 $t\to\infty$ 均为零，确得原接触常数。于是离散余弦和除以 $2N$ 恰是中点多项式减极化多项式，半径为 $e^{-x/N},e^{-y/N}$；再用已证同半径 Hermite 下界。若要得到无未证前提的 Lean 声明，必须由真实连续谱的余量定理驱动这个离散步骤，不能直接引用 Python 输出。
