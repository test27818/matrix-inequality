# 原核在浅端紧方块的有限 κ 二阶差商重建

**范围。** 对 $0\le\kappa\le1/16$、$0\le x\le y\le3$，原频谱在任意实频率的中点 Hermite 与极化对数核之差的余弦变换非负。特别地，这覆盖所有整数 $n\ge16$ 的原离散共用核在此方块的全部角度。并非整个 $[0,32]^2$；也不以早期 C++ 留言或日志作证据。

## 降维公式

用 [原谱定义、归一化与离散取样](moderate-helmholtz.md)。记 $q=1-p$、$b=\kappa^2$、$P=q(q^2-b)$、$z=q^2$、$C(v)=(\cosh v-1)/v^2$。注意 $q$ 不是 $\kappa$。在 $\kappa<q<1$ 有
$$
 F_p(v)=P f(q,v),\qquad
 f(q,v)=\frac v{C(v)}\sum_{a,h\ge0}\frac{(qv)^{2a}(\kappa v)^{2h}}{(2a+2h+3)!}.
$$
这是两个双曲正弦级数的差商，对 $q=\kappa$ 仍有上述无奇点延拓。令 $U_t=\int_0^t f(q,v)\,dv$，$V=\int_0^\infty f(q,v)\,dv$。展开 $(\cosh v-1)^{-1}=2\sum_{j\ge1}j e^{-jv}$，对 $G_a(v)=\sinh(av)/a$ 逐项积分；$\int_0^\infty e^{-jv}G_a(v)\,dv=(j^2-a^2)^{-1}$，故
$$
 \Phi=(1-z)V=\frac2{1-b}
   +2\sum_{j\ge2}\frac{j^2}{j^2-b}\frac{1-z}{j(j^2-z)}.
$$
置 $K_t=\Phi-(1-z)U_t$、$m=(x+y)/2$、$d=y-x$、$v_u=x+du$、$\lambda(u)=\min(u,1-u)$，
$$
 W=\int_0^1 f(q,v_u)\,du,\quad
 Y=\frac{W^2}{2(\sqrt{K_x}+\sqrt{K_y})^2},\quad
 B=\frac12\int_0^1\lambda(u)\partial_v f(q,v_u)\,du.
$$
$K_t$ 是正数在参数端点处的连续延拓；下述程序同时检查所有被使用的平方根和分母均为正。中点 Peano 核作用于 $U_t$ 给
$(U_x+U_y)/2-U_m=d^2B$。再用 $(a+b)/2-\sqrt{ab}=(a-b)^2/[2(\sqrt a+\sqrt b)^2]$：
$$
 D(p)=J_m(p)-\sqrt{J_x(p)J_y(p)}
   =d^2P\{B+(1-z)Y\}. \tag{1}
$$
右边包含 $d=0$ 与 $p=0$，均无相近半径相减或 $1/p$ 的不定式。

**与连续极限共用的主项。** 因为原 $F_p(v)$ 分子中的 $G_\kappa$ 项为 $q$ 的一次式，
$\partial_q^2 F_p(v)=v^2\sinh(qv)/(\cosh v-1)$ **与 $b$ 无关**。令
$S(v)=vC'(v)$、$\operatorname{sinhc}(t)=\sinh(t)/t$，
$$
 H(q,v)=\frac{\cosh(qv)-(S(v)/C(v))\operatorname{sinhc}(qv)}{C(v)}.
$$
直接微分 $v^2\sinh(qv)/(q(\cosh v-1))$ 得 $H$，故 $P_0=(PB)''/q=\frac12\int\lambda H$。所有撇號在本节为 $q$ 微分。于是
$$
 \frac{D''(p)}{q d^2}
 =P_0 +[6(1+b)-20z]Y+
 [6(1+b)q-10qz-2b/q]Y'+[(1+b)z-z^2-b]Y''. \tag{2}
$$
(2) 在 $q=0$ 按连续极限取值；在物理区 $0\le b/q\le\min(q,1/16)$，检查器不对含零区间做除法。

**有限 $b$ 的控制。** 设 $E(t)=\sum_{a\ge0}t^{2a}/(2a+3)!$；对 $j=0,1,2$，$\partial_q^j f$ 均非负且夹在 $v^{j+1}E^{(j)}(qv)/C(v)$ 和其乘上 $1/(1-bv^2/20)$ 之间。这来自逐系数比
$(2a+3)!/(2a+2h+3)!\le20^{-h}$，同样适用于正项微分，故不是把 $b=0$ 的数值替代有限 $b$。特别地由 $3(\sinh v-v)\le v(\cosh v-1)$ 得 $f\le(v/3)/(1-bv^2/20)$。$\Phi,\Phi_z,\Phi_{zz}$ 的尾从 $j=N+1$ 开始，绝对值分别不超过
$g/[N(N+1)]$、$g/[N(N+1)]$、$2g/[N(N+1)((N+1)^2-1)]$，其中 $g=(N+1)^2/((N+1)^2-b)$，由望远镜求和 $\sum_{j>N}1/[j(j^2-1)]=1/[2N(N+1)]$ 即得。检查器取 $N=48$。

**有限截止原子。** 不能误用 $\kappa=0$ 的端点 $D'(1)=0$。(1) 给
$$
 D'_p((1-\kappa)-)=-2b d^2\{B+(1-b)Y\}_{q=\kappa}.
$$
在端点的一组独立闭盒上验证大括号非负。记 $E_* =\sum_{a,h}(qv)^{2a}(\kappa v)^{2h}/(2a+2h+3)!$；因为 $\partial_vE_*\ge0$，
$\partial_v f\ge (1-S/C)E_*/C$。以 $q\in[0,1/16]$ 的 $\Phi$ 上界 $K_{\max}$ 同时控制 $K_x,K_y$，则
$B+(1-b)Y\ge\tfrac12\int\lambda\partial_v f+(1-b)W^2/(8K_{\max})$。沿 $u$ 的 16 个等长子区间用区间下端积分这个表达式，负号因子使用正确的外包乘法。结合 $D''\ge0$、$D(1-\kappa)=0$ 与端点非正斜率，$D$ 连续延拓为截止后的零函数，二阶导数为非负测度。在 $\omega\ne0$ 时，余弦积分等于 $\omega^{-2}\int(1-\cos(\omega p))\,dD'(p)\ge0$；$\omega=0$ 时由凸性与端点得 $D\ge0$。$d=0$ 时 $D\equiv0$。由中深窗文档第5节所给取样论证得到所有离散角度的原核比较。

## 精确算术与重放

[实现](../../certificates/verify_finite_shallow.py) 只复用 [连续极限版](../../certificates/verify_continuum_shallow.py) 中的整数有界算术、正系数 Taylor 界、平方根 jet。变量空间的正方格、三维自适应细分、径向加权积分写在新脚本中。以 $2^{56}$ 为分母作外舍入，$e$、$\pi$ 的浮点实现均未参与。所有正系数函数的 48 阶 Taylor 余项小于 $2\cdot4^{49}/49!<2^{-56}$（由逐项比控制）。在 $u$ 分块上对 $W,W',W''$ 作区间，两个端点间的全径向格可用正项前缀和收紧；对 $\lambda$ 用积分的精确有理权重。所有在正方格中保留的交集都是同一平均数的多个外包。只对下端尚非正的盒作三参数二分；无序半径由对称性交换回 $x\le y$。
 
    PYTHONDONTWRITEBYTECODE=1 python3 certificates/verify_finite_shallow.py 3 8 16 1 7

[整数日志：从源码重放](../../replay.sh)给出每层搜索尺寸及最后严格余额。结论为
$D''(p)\ge (1157372000/2^{56})(1-p)(y-x)^2$，截止大括号下端为 $393048156776508/2^{56}$。这里是解析证明加可重复运行的有理区间证书，不是全问题的 Lean 闭合。
