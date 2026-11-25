# 浅－深有限交叉窗的极点消去检查器

独立重建三个历史闭窗 $x\in[3,9],[9,12],[12,16]$、$y\in[32,48]$、$b=\kappa^2\in[0,1/256]$。这里的数学降维与 [中深交叉窗](moderate-helmholtz.md) 共用真实 $J_t,T_t,W_t,D$，但更换分点并显式保留浅端前 3 个 alias。旧 `psd_moderate_48_cross.cpp` 仍未恢复；这不是旧代码或旧成功日志的复制，也不是本区域的 Lean 重放。

## 高频段

取 $\beta=3/4<1-\kappa$，$c=(y-x)/2$，$m=(x+y)/2$。正卷积支撑的单调性、凸性及截止原子吸收不变。径向积分的 $v\ge x\ge3$，故 $p v\ge\beta\cdot3>1$，仍有 $0\le F_p(v)\le e^{-pv}$。所以高频项和截止原子整个下界仍为 $2V'(\beta)-c^2e^{-\beta m}/(\beta m)$，其中 $V=\sqrt{J_x J_y}$。不对截止前的连续高频段做采点。

## 低频段的浅端余尾

设 $q=e^{-t}$，对 $t\ge10$ 用 [原中深版](moderate-helmholtz.md) 的一个 alias。若某参数盒可含 $t<10$，取 $K=3$ 并显式保留深度 $t$ 的前 $K$ 个 alias；这一选择只决定 enclosure 算法，不改变目标函数。

以 $M_i(t,a)=\int_t^\infty v^i e^{-av}/i!\,dv$ 表示单个指数尾矩。在 $e^t T_t$ 的两个移动指数列对第 $i=0,1,2$ 槽保留
$$
 e^t\sum_{j=1}^K j M_i(t,j+1-p),\qquad e^t\sum_{j=1}^K (j+1)M_i(t,j+p).
$$
对于 $v\ge t$，相对首 alias 的余下权重和逐点不超过
$$
 \sum_{r=K}^\infty(r+2)q^r=\frac{q^K[(K+2)-(K+1)q]}{(1-q)^2}
 \le \frac{(K+2)q^K}{(1-q)^2}.
$$
因此乘以相应首个 $e^t M_i$ 就给全正余尾上界。先按正尾估计，最后按 $T,T',T''/2$ 原来符号加减。

常数双曲列 $e^t H_t$ 各 alias（含权重）保留正 Taylor 展开的 $j=0,\ldots,23$：
$$
 2l\,e^{(1-l)t}\sum_{j=0}^{23} b^j F_{2j+1}(t,l),
 \qquad F_n(t,l)=\frac{1}{l^{n+1}}\sum_{i=0}^n\frac{(lt)^i}{i!}.
$$
实现递推 $F_0=1/l,\ F_{n}= (F_{n-1}+t^n/n!)/l$。单个加权 alias 被删的 $b$ 尾严格小于 $2^{-100}$，因为 $e^t M_{2j+1}(t,l)\le e^t/l^{2j+2}$，$l\ge1$ 且沿用 $t\le56$ 的旧整数估计。未保留的所有 $l>K$ 双曲 alias 相对完整首 alias 不超过 $(K+1)q^K/(1-q)^2$。该上界和已付的每个 $b$ 尾叠加不重算遗漏项。

## 复跑

```sh
PYTHONDONTWRITEBYTECODE=1 python3 certificates/shallow_cross_helmholtz.py --root all --mesh 6 \
  --export certificates/shallow-cross-helmholtz-witness.json
```

已两次独立复跑通过：[64 单元：从源码重放](../../replay.sh)（并导出[精确参数分割：从源码重放](../../replay.sh)）和 [128 单元：从源码重放](../../replay.sh)。两次参数树计数一致。

| $x$ 范围，均配 $32\le y\le48$ | 访问 | 叶 | 最大深度 | 64 单元最小余量分子，公分母 $2^{112}$ |
|:--|--:|--:|--:|--:|
| $[3,9]$ | 155 | 78 | 9 | 127357935259487969634053478486322782 |
| $[9,12]$ | 1 | 1 | 0 | 108040726599784003621575550333505474 |
| $[12,16]$ | 3 | 2 | 1 | 40184919039754252119241389644139768 |

整数端点精度、外舍入、jet 广义四则及参数树遍历与中深版一致。统一谱格每个端点确取 $3j/(4K_p)$，$K_p=2^{\mathrm{mesh}}$。三个根都对整个 $b$ 区间接受，且数值都对应正因子 $e^x$ 缩放后的连续余弦余量，最小严格大于 7。将这些根送入现有 Lean 重放仍须换用 $\beta=3/4$ 和三-alias 版本。ShallowCrossAlias.lean 已形式化一般 $K$ 的正移动指数尾界，未把 JSON 的标签当成 Lean 定理。

**单-alias 备选。** 加 `--single-alias` 已用同一 $\beta=3/4$ 跑通全部根（访问/叶 $627/314,1/1,3/2$）。[日志：从源码重放](../../replay.sh)、[精确分割：从源码重放](../../replay.sh) 分别保留。因此 Lean 的另一路径只需更换固定谱分点，不必重写实 alias 语义；ShallowCrossEval.lean 为这个参数版本的可执行定义，尚未声称其整条真实谱桥梁已过内核。
