# 中小深度紧窗：连续核的双参数图与全频率验证

这里证明的是 **$\kappa=0$ 的原中点 Hermite 混合差**，不是 $\kappa>0$。覆盖整块参数的结论仅以本文列出的 `PASS` 日志为准。与原矩阵的连接使用根提升后取连续极限，不用把这里的严格符号冒认为每个有限阶原证书的符号。

取 $p\in[0,1]$，$q=1-p$，$J_v=\int_v^\infty(\sinh(qs)-qs)/(\cosh s-1)\,ds$，$m=(x+y)/2,h=(y-x)/2$。证明 $D(p)=J_m-\sqrt{J_xJ_y}$ 在 $p>1$ 接零后的余弦变换全频非负。归一化目标 $d=e^m D/h^2$，$h=0$ 取极限（原差恰为零）。两种参数图均给出同一 $d$ 及 $\partial_q^2d=d''(p)$ 的包络。

## 1. 在靠近 $p=0$ 处不作极点相减

置 $B(v)=(1-p e^{pv}J_v)/p$，$W(v)=1-pB(v)$，$Y_r(v)=e^{qv}\partial_v^r B(v)$。用 $(\cosh v-1)^{-1}=2\sum_{j\ge1}je^{-jv}$ 逐项积分，并先扣去 $1/p$ 后，得无奇点的式子。设 $a=j-1+q,k=j-1+2q$，对 $r=0,1,2$，
$$
Y_r=\sum_{j\ge1}e^{-(j-1)v}\left\{
2q\left[(-a)^r(v+1/j)+r(-a)^{r-1}\right]
 +\frac{j(-k)^r}{j+q}e^{-qv}-\frac{(j+1)(-j)^r}{j+1-q}e^{-(1-q)v}\right\}.
$$
$r=0$ 时方括号的第二项定义为零。令 $v_s=m+hs$，
$$
 U=\frac12\int_{-1}^1e^{-qhs}Y_1(v_s)\,ds,\qquad
 V=\frac12\int_{-1}^1(1-|s|)e^{-qhs}Y_2(v_s)\,ds.
$$
即 $U=e^{qm}(B(y)-B(x))/(2h)$、$V=e^{qm}((B(x)+B(y))/2-B(m))/h^2$；径向相近时由积分定义，不相减。把算术与几何均值之差写成平方差除法得到
$$
 d= V+2p e^{-qm}\frac{U^2}{(\sqrt{W(x)}+\sqrt{W(y)})^2},
 \qquad W(v)=1-p e^{-qv}Y_0(v).
$$
此式在 $p=0$ 仍无奇点。在 $h$ 离零时，区间也可直接用上述端点商；两个界同时计算时取交集。

## 2. 在 $p=1$ 附近先提出精确三次零

令 $A(v)=e^vJ_v/q^3$，在 $q=0$ 取连续极限。若 $E_k(t)=\sum_{i=0}^kt^i/i!$，$E_{-1}=0$，$C_k(t)=\sum_{i\ge0}t^{2i}/(2i+k)!$，则
$$
 P_d(j)=\sum_{i\ge0}(q/j)^{2i}E_{2i+d}(jv)
 =\frac{E_{d-2}(jv)+(jv)^{d-1}C_{d-1}(qv)+(jv)^d C_d(qv)}{1-(q/j)^2}.
$$
由对 $\sinh(qv)-qv$ 的正项展开和指数矩的积分，$r=0,1,2$ 有
$$
 A^{(r)}(v)=2\sum_{j\ge1}j^{-3}e^{-(j-1)v}
 \sum_{t=0}^r {r\choose t}(1-j)^{r-t}j^t P_{3-t}(j).
$$
只在 $q\le1/2$ 使用，故分母远离零。定义
$$
U=\tfrac12\int_{-1}^1 A'(v_s)\,ds,qquad
V=\tfrac12\int_{-1}^1(1-|s|)A''(v_s)\,ds.
$$
同一算术-几何均值恒等式给
$$
 d=q^3\left(2U^2/(\sqrt{A(x)}+\sqrt{A(y)})^2-V\right).
$$
这同时证明 $d(1)=d'(1)=0$。在中间可使用两图任一实际通过的包络；相互重叠只提高精度。

## 3. 绝对收敛、级数余项和舍入

在 $3\le v\le24$ 取最小正整数 $n$ 使 $n v_{\min}\ge105$，保留 $j=1,\ldots,n$；$n\le35$。在复圆盘 $|q'-q|\le1/8$ 中，第一图中的每项（除 $e^{-(j-1)v}$）在 $r\le2$ 时模不超过 $1000(j+2)^2$：有 $|q'|,|1-q'|\le9/8$，移动指数模不超过 $e^3<21$，两个有理系数分别不超过 $8/7,16/7$。第二图中 $|q'|\le5/8$，$e^{|q'|v}<4\cdot10^6$，$|P_d|<10^{12}j^3$；连 $2j^{-3}$ 与 $r$ 次径向导数后不超过 $8\cdot10^{12}j^2$。Cauchy 不等式对 $q$ 的前三个 Taylor 槽（实际二阶导数的一半是第三槽）最多再乘 $64$。因此共同的保守尾绝对值
$$
 10^{16}\sum_{j>n}(j+2)^2e^{-(j-1)v}\le 2\,10^{16}38^2 e^{-105}<2^{-45}.
$$
$C_k$ 保留 $i=0,\ldots,34$，全用正项；因 $qv\le12$，其三个 Taylor 槽的共同余项至多 $2\cdot24^2\,12^{68}/68!<2^{-45}$。大常数的比较用 [有理数短脚本](../../certificates/compact_continuum_mixed_bounds.py)独立验证。

程序用 IEC-559 binary64，每次加、减、乘、除均把候选值上下各拓一格，不倚赖未给误差界的库指数。非负指数自变量减半至 $\le1/16$，求 12 次 Taylor 和，余项由两倍第13项向上包络，再反复平方；负指数取倒数。平方根候选用向外取整的平方进行回证。程序要求默认就近舍入，禁止 fastmath 与融合乘加。径向区间在每个小格中覆盖端点，合并若干格用上下确界的外包；$s$ 积分取均匀子格，对 $1-|s|$ 用准确的三角权重；直接端点差与 Peano 表达互不相容必须中止，而不会当作通过。径向平方根无法验证正分母时，该参数图不能用于接受该单元。

## 4. 全实频率，不是离散测试

每个 $p$ 格得到 $d,d''$ 的实区间，最终以 $2^{-28}$ 为量子向外取整成整数。首先可用[带符号曲率下包络](signed-curvature.md)（此处截止原子零）。另有更紧的双侧判定。把格上曲率区间取一个整数中心 $c_i$ 和误差半径，值区间取中心 $f_i$ 和误差半径；$E,F$ 分别记格宽加权的半径和。令 $\delta=-d'(0)$。用 [去极点端点展开](../continuum/diagonal.md)得
$$
\delta=\frac{e^m}{2}\bigl(1-[\operatorname{avg}_{[x,y]}C'(v)]^2\bigr)
 =\frac12 (2-\operatorname{avg}a)(e^m\operatorname{avg}a),\quad
 C(v)=v\coth(v/2),\quad a=1-C'(v)=2e^{-v}(v-1+e^{-v})/(1-e^{-v})^2.
$$
后者用 16 个闭盒分别包络两个平均，不减近数。当 $h=0$ 即单点评价。令 $V=|c_{N-1}|+\sum_{i>0}|c_i-c_{i-1}|$（按 $p$ 递增）。对 $\omega>0$ 两次分部积分给 $\omega^2\widehat d=\delta-\int d''\cos$；按格积分再差分给 $\omega^2\widehat d\ge\delta-E-V/\omega$，因此严格护栏以上的整条无界频率均被覆盖。

所有物理参数上的 $d\ge0$ 还由[径向严格对数凹性的正系数论证](../continuum/diagonal.md#proven-mixed-low-frequency-strip)给出，故对值区间可把负下端提高到零。较小频率用 $0\le\omega\le1$ 上 $1/2\le\cos(\omega p)\le1$；其后把频率每 $1/2$ 一格。对每个频率格，用[同一精确整数余弦权重](../../certificates/signed_curvature_fast.cpp)包围 $\alpha_i$，即 $\cos(\omega p)$ 在第 $i$ 个 $p$ 格的平均。两个量
$$
 N^{-1}\sum_i f_i\alpha_i-F,\qquad
 \delta-N^{-1}\sum_i c_i\alpha_i-E
$$
分别是 $\widehat d$、$\omega^2\widehat d$ 的下界；任一个在整频率格都严格为正即可。高频阈值超出上限或某格失败只会要求径向二分，不表示反例。对径向子盒不可按浮点中点接受；此实现使用 dyadic 径向与谱网格，参数盒验不出则对两径向各二分，已接受的盒免于重算。

## 5. 重放

以下是 κ=0 中点差，$x\le y$，$x=y$ 原差恰等于零。

```sh
cd certificates
python3 compact_continuum_mixed_bounds.py
g++ -std=c++20 -O3 -ffp-contract=off -fno-fast-math -fsanitize=undefined -fno-sanitize-recover=all compact_continuum_mixed.cpp -o /tmp/compact-replay
/tmp/compact-replay 8 12 8 16 2 64 5 8
/tmp/compact-replay 3 8 3 16 2 64 7 16
```
已通过的区域合起来是 **$3\le x\le12,\ x\le y\le16$**。$8\le x\le12$ 见[中间闭域重放：从源码重放](../../replay.sh)，$3\le x\le8$ 见[浅侧闭域重放：从源码重放](../../replay.sh)：浅段五层分别解出 0、91、2074、2682、4424 个径向盒，最后一层没有未决格。两段均未使用同半径的额外对角松弛。命令尾四项为每单位初始径向格个数、初始频谱格数、最多层数、深层 Peano 格数（浅层依半间距可少用格）。失败的可选点或未完成的输出不是区域证明。

两条命令的最新独立 UBSan 实跑均以退出码 0 结束，见[中段独立日志：从源码重放](../../replay.sh)和[浅段独立日志：从源码重放](../../replay.sh)。一键命令包含纯有理常数复核、严格编译模式，以及两个主根的完整参数。独立数学审阅逐项核对了两张参数图、相近半径的 Peano 因子、复 $q$ 的尾界、截止处的零斜率、整条频率轴的保证以及 $q,p$ 网格倒序。
