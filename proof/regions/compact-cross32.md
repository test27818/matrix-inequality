# 原 [0,32] 方窗的三条离对角上沿闭盒

采用 [浅深交叉窗](shallow-cross-helmholtz.md) 已推导的原 $F,J,D$、48 项双曲 Taylor 函数、高频截止原子吸收、正余尾及 [中深 Helmholtz](moderate-helmholtz.md) 公式 (1)–(6)，仅放宽参数盒为 $3\le x\le16$, $x<y\le32$。这些推导的实际边界依赖是：
* $v\ge t\ge3$ 与 $\beta\ge1/2$，故 $2e^{-\beta v}\le2e^{-3/2}<1$；高频界完整有效；
* $\beta\le29/32<1-1/16$，截止始终在严格高频段内部；
* $t\le32<56$，指数矩的 Taylor 余界并不增加；其余收敛与正性只用 $t>0$、$b\le1/256$；
* $\beta/2^m$ 在本重放均为能够精确表示的 dyadic，完整闭网格不漏区间。

因此 [复用整数实现](../../certificates/compact_cross_helmholtz.py) 的每个接受闭叶均给 $e^x(c^2+\omega^2)\int_0^{1-\kappa}D(p)\cos(\omega p)\,dp>0$，不是只在网格端点验算。最后如中深文档第5节离散取样。$b=\kappa^2\in[0,1/256]$，三块均不限制角度。

    python3 certificates/compact_cross_helmholtz.py 3 9 24 32 3/4 6
    python3 certificates/compact_cross_helmholtz.py 9 12 24 32 7/8 6
    python3 certificates/compact_cross_helmholtz.py 12 16 24 32 29/32 7

[整数日志：从源码重放](../../replay.sh)分别得到 40、16、45 个通过的闭叶，严格下界分母均为 $2^{112}$。未由这些程序覆盖的相邻近对角盒不能引用本页。  
