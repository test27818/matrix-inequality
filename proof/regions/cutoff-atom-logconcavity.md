# 截止点原子的全局符号

设 $0<\kappa\le1/4$，$q=1-p$，$L=1-\kappa$。取 $P=q(q^2-\kappa^2)$、
$A_z=J_z/P$。根据[正卷积公式](moderate-helmholtz.md)，$A_z$ 在 $q=\kappa$ 有连续可微的正延拓。
则对于**所有**深度 $0\le x\le y$、中点 $m=(x+y)/2$，
$$D'(L-)=-2\kappa^2(A_m-\sqrt{A_xA_y})\big|_{q=\kappa}\le0.$$
所以分布曲率 $D''$ 在右端的原子 $-D'(L-)$ 非负，在 Helmholtz 符号测度的预算中不必再减去正斜率。$\kappa=0$ 时 $P=q^3$ 给零原子；两个结论都不假设中点缺陷的余弦变换已证。

**证明。** $G_q(v)=\sinh(qv)/q$，取关于 $q^2$ 的导数可得 $A_z(L)$ 是 $v\in[z,\infty)$ 上正密度
$$
\frac{h(\kappa v)}{2\kappa^3(\cosh v-1)},\qquad h(t)=t\cosh t-\sinh t
$$
的尾积分。这个密度关于 $v>0$ 严格对数凹。事实上
$$
(h')^2-hh''=\sinh^2t-t^2\ge t^4/3,
\qquad 0<h(t)=\int_0^t u\sinh u\,du\le (t^3/3)\cosh t,
$$
第一式用 $\sinh t\ge t+t^3/6$，第二式用 $\sinh u\le u\cosh t$。因此
$$
-\partial_v^2\log h(\kappa v)\ge\frac{3}{v^2\cosh^2(\kappa v)}.
$$
另一方面 $\sinh(v/2)=2\sinh(v/4)\cosh(v/4)\ge(v/2)\cosh(v/4)$，所以
$$
-\partial_v^2\log(\cosh v-1)=\frac{1}{\cosh v-1}
\le\frac{2}{v^2\cosh^2(v/4)}\le\frac{2}{v^2\cosh^2(\kappa v)}.
$$
两式相减即所需严格负二阶导。正对数凹密度 $f$ 的右尾 $A_z$ 也对数凹：对每个 $s>0$ 比值 $f(z+s)/f(z)$ 非增，从而 $A_z/f(z)=\int_0^\infty f(z+s)/f(z)\,ds$ 非增、$(\log A_z)'=-f(z)/A_z$ 亦非增。对 $z=0$ 取右极限。于是 $A_m\ge\sqrt{A_xA_y}$。最后对 $D=P(A_m-\sqrt{A_xA_y})$ 从内段向 $p=L$ 求导，$P(L)=0$, $P'(L)=-2\kappa^2$，得证。

注意：本引理证明的是**一个原子**的符号，不自动证明内段曲率、内段中点缺陷之号，亦不自动证明全频率余弦变换之号。
