# PSD Matrix Product Inequality | 半正定矩阵乘积不等式

A reproducible proof of a product inequality for **complex Hermitian positive semidefinite (PSD) matrices**, using Fourier analysis and rigorous interval arithmetic.

**半正定矩阵不等式。** 对任意正整数 $n$ 与复 Hermitian 半正定矩阵 $A=(a_{ij})\in\mathbb C^{n\times n}$，都有

$$
\prod_{i,j=1}^{n}\lvert 1-a_{ij}\rvert
\;\le\;
\prod_{i=1}^{n}\lvert 1-a_{ii}^{\,n}\rvert.
$$

The product on the left uses all ordered pairs, including the diagonal.

## Start here | 阅读入口

- **[完整证明 / Full proof (中文)](proof/PROOF.md)**：从矩阵归约到连续标量核，再到各参数区间的完整覆盖。The main write-up is in Chinese and links to the relevant analytic lemmas and certificate specifications.
- **[Analytic lemmas / 分区引理](proof/)**：按连续归约与参数区域组织。
- **[Verification source / 证书源码](certificates/)** 和 **[replay.sh](replay.sh)**：重建并运行证明所用的核验。

**Proof status / 证明状态。** The complete proof is computer-assisted. The reduction uses Schur products and Fourier positivity; rigorous outward-rounded interval certificates cover the remaining scalar inequalities. These checks bound whole parameter intervals, not just sampled points. The public proof does not require Lean.  
完整结论已证。区间证书是证明的一部分，不是数值实验；本仓库不需要 Lean。

## Reproduce the checks | 重放核验

从仓库根目录运行：

~~~sh
sha256sum -c SHA256SUMS
sh replay.sh all
~~~

Requires Python 3, a C++20 compiler, GMP development headers, the MPFR 6 shared library, and an IEEE 754 binary64 environment. 依赖 Python 3、C++20 编译器、GMP 开发头文件、MPFR 6 共享库及 IEEE 754 binary64 环境。详细分区见[证明正文](proof/PROOF.md)。

If you cite the proof or report a verification issue, include the commit SHA so readers can reproduce the exact version.
