# 半正定 Hermitian 矩阵的乘积不等式

对任意 $n$ 阶复 Hermitian 半正定矩阵 $A=(a_{ij})$，有
$$
\prod_{i,j=1}^{n}|1-a_{ij}|\le\prod_{i=1}^{n}|1-a_{ii}^{n}|.
$$

- [证明正文](proof/PROOF.md)：矩阵归约、连续标量核、完整覆盖。
- [proof/](proof/)：按用途和参数区间组织的解析引理。
- [certificates/](certificates/)：重放数学证书所需的源码及最小数据。
- [replay.sh](replay.sh)：从源码重建并重放全部证书。

在本目录运行 `sh replay.sh all`。需要 Python 3、C++20、GMP 开发头文件、
GNU MPFR 6 共享库以及 IEEE binary64 算术。区间核验为计算机辅助数学证明
的一部分，读者可从源码独立重放；不依赖归档日志。
