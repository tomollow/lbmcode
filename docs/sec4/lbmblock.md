# lbmblock.c 説明ドキュメント

## 概要

[src/sec4/lbmblock.c](../../src/sec4/lbmblock.c) は、2 次元 Couette 流れを **多重格子（Multi-block）格子ボルツマン法** で解くコードです。領域を上下に分け、上半分を **粗格子**（coarse grid）、下半分を **細格子**（fine grid, 格子幅 1/2）で覆い、両格子をインターフェースで接続して時間進行させます。

主な技術要素：

- **D2Q9 + BGK 衝突**（[lbmcm.md](lbmcm.md) と同じ平衡分布）
- **格子間で緩和時間を再スケーリング**して物理動粘性 $\nu$ を一致させる
- **インターフェース転送**で粗・細格子の境界分布関数を交換：
  - **細 → 粗**: 非平衡部のスケーリング（粒子分布の再構築）
  - **粗 → 細**:
    - **Lagrangian 2 次時間補間**（粗が 2 細ステップに 1 回しか進まないため、$n-\tfrac{1}{2}$ 時刻を内挿）
    - **3 次スプライン空間補間**（細の x 解像度が 2 倍のため、中間点を埋める）
- **Zou-He 非平衡 bounce-back** で上下壁条件
- **x 周期境界**

検証対象は線形 Couette 解析解 $u(y) = u_t \, y / H$ で、コードは時間 100, 200, 500, 1000, 3000 の中心断面 $u$ プロファイルをファイル出力します。

## 多重格子の配置

格子定数（[lbmblock.c:90, 108-111](../../src/sec4/lbmblock.c#L90-L111)）：

```c
nx = 20, ny = 32, m = 2;
nxf = m*nx = 40;  nyf = ny = 32;       // fine grid
nxc =   nx = 20;  nyc = ny/m + 1 = 17; // coarse grid
```

物理座標 $y_{\rm phys} \in [0, n_y]$ における各格子のカバー範囲（refinement ratio $m = 2$）：

| 領域 | 物理 $y$ 範囲 | 格子 | 格子幅 |
|---|---|---|---|
| 下半分（壁面側） | $y_{\rm phys} \in [0, 16]$ | 細格子 $y_f \in [0, 32]$ | $\delta x_f = 1/2$ |
| 重なり（オーバーラップ） | $y_{\rm phys} \in [15, 16]$ | 両方 | — |
| 上半分（蓋側） | $y_{\rm phys} \in [15, 32]$ | 粗格子 $y_c \in [0, 17]$ | $\delta x_c = 1$ |

インターフェースの対応関係：

- **粗 $y_c = 0$** $\Leftrightarrow$ **細 $y_f = n_{y,f} - 2 = 30$**（共通の物理位置 $y_{\rm phys} = 15$）
- **粗 $y_c = 1$** $\Leftrightarrow$ **細 $y_f = n_{y,f} = 32$**（共通の物理位置 $y_{\rm phys} = 16$）
- **粗 $y_c = n_{y,c} = 17$**: 上端の動く壁 $u_t$
- **細 $y_f = 0$**: 下端の固定壁 $u_b$

x 方向は両格子ともに周期境界、ただし細格子は粗格子の 2 倍の解像度（$2i_c \leftrightarrow i_f$）。

## Couette 流れと解析解

定常 Couette 流れの解析解（[lbmblock.c:562-566](../../src/sec4/lbmblock.c#L562-L566) の参照値）：

$$
u(y) = u_b + \frac{u_t - u_b}{H}\,y
$$

本実装では $u_b = 0$, $u_t = 0.01$, $H = n_y = 32$ なので

$$
u(y) = u_t\,\frac{y}{n_y}
$$

最終ステップで標準出力に格子中心断面 $u(n_x/2, j)$ と解析解の比較が表示されます。

## D2Q9 と BGK 衝突

[lbmcm.md](lbmcm.md) と同じ D2Q9 格子（[lbmblock.c:122-126](../../src/sec4/lbmblock.c#L122-L126)）。平衡分布 $f_k^{\rm eq}$ も同形：

$$
f_k^{\rm eq} = w_k\,\rho\left[1 + 3\,(\mathbf{c}_k\!\cdot\!\mathbf{u}) + \tfrac{9}{2}(\mathbf{c}_k\!\cdot\!\mathbf{u})^2 - \tfrac{3}{2}|\mathbf{u}|^2\right]
$$

粗・細格子の各セルで BGK 衝突（[lbmblock.c:216-218, 315-317](../../src/sec4/lbmblock.c#L216-L218)）：

$$
f_k^{\rm post,c} = f_k^c - \frac{1}{\tau_c}(f_k^c - f_k^{\rm eq,c}),\qquad
f_k^{\rm post,f} = f_k^f - \frac{1}{\tau_f}(f_k^f - f_k^{\rm eq,f})
$$

## 緩和時間スケーリング（鍵となる関係）

LBM の動粘性 $\nu$ は格子単位で

$$
\nu = c_s^2\left(\tau - \tfrac{1}{2}\right)\frac{\delta x^2}{\delta t} = \tfrac{1}{3}\left(\tau - \tfrac{1}{2}\right)\frac{\delta x^2}{\delta t}
$$

均等細分（$\delta x_f = \delta x_c / m$, $\delta t_f = \delta t_c / m$）で **物理動粘性を両格子で一致させる** には：

$$
\nu_{\rm phys}
= \tfrac{1}{3}(\tau_c - \tfrac{1}{2})\frac{\delta x_c^2}{\delta t_c}
= \tfrac{1}{3}(\tau_f - \tfrac{1}{2})\frac{\delta x_f^2}{\delta t_f}
= \tfrac{1}{3}(\tau_f - \tfrac{1}{2})\frac{\delta x_c^2}{m\,\delta t_c}
$$

これより

$$
\boxed{\;\tau_f = \tfrac{1}{2} + m\left(\tau_c - \tfrac{1}{2}\right)\;}
$$

コードでは [lbmblock.c:111](../../src/sec4/lbmblock.c#L111)：

```c
tauf = 0.5 + (double)m*(tauc - 0.5);
```

$\tau_c = 1.2$, $m = 2$ で $\tau_f = 0.5 + 2 \times 0.7 = 1.9$。

格子単位での動粘性は [lbmblock.c:114-115](../../src/sec4/lbmblock.c#L114-L115)：

$$
\nu_c = \frac{\tau_c - 1/2}{3} = 0.2333,\qquad
\nu_f = \frac{\tau_f - 1/2}{3\,m} = 0.2333
$$

両者が一致することで、Chapman-Enskog 展開で得られる NS 方程式の粘性係数が物理的に連続になります。

## 壁面境界条件（Zou-He 非平衡 bounce-back）

### 上端の動く壁（粗格子 $y_c = n_{y,c}$, [lbmblock.c:271-279](../../src/sec4/lbmblock.c#L271-L279)）

壁面の密度を既知の分布関数から評価：

$$
\rho_{\rm top} = \frac{f_0 + f_1 + f_3 + 2(f_2 + f_5 + f_6)}{1 + v_t}
$$

未知の壁向き分布 $f_4, f_7, f_8$ を Zou-He 法で構築：

$$
\begin{aligned}
f_4 &= f_2 - \tfrac{2}{3}\rho_{\rm top}\,v_t \\
f_7 &= f_5 + \tfrac{1}{2}(f_1 - f_3) - \rho_{\rm top}\!\left(\tfrac{u_t}{2} + \tfrac{v_t}{6}\right) \\
f_8 &= f_6 - \tfrac{1}{2}(f_1 - f_3) - \rho_{\rm top}\!\left(-\tfrac{u_t}{2} + \tfrac{v_t}{6}\right)
\end{aligned}
$$

本実装では $u_t = 0.01$, $v_t = 0$。

### 下端の固定壁（細格子 $y_f = 0$, [lbmblock.c:369-377](../../src/sec4/lbmblock.c#L369-L377)）

対称的に、未知の $f_2, f_5, f_6$ を上向き運動量を考慮して構築：

$$
\rho_{\rm bot} = \frac{f_0 + f_1 + f_3 + 2(f_4 + f_7 + f_8)}{1 - v_b}
$$

$$
\begin{aligned}
f_2 &= f_4 + \tfrac{2}{3}\rho_{\rm bot}\,v_b \\
f_5 &= f_7 - \tfrac{1}{2}(f_1 - f_3) + \rho_{\rm bot}\!\left(\tfrac{u_b}{2} + \tfrac{v_b}{6}\right) \\
f_6 &= f_8 + \tfrac{1}{2}(f_1 - f_3) + \rho_{\rm bot}\!\left(-\tfrac{u_b}{2} + \tfrac{v_b}{6}\right)
\end{aligned}
$$

本実装では $u_b = v_b = 0$。

## インターフェース転送

### 非平衡部のスケーリング則

両格子で **平衡分布** は連続だが、**非平衡分布** $f_k^{\rm neq} = f_k - f_k^{\rm eq}$ は緩和時間が異なるため不連続。Chapman-Enskog 展開から、応力テンソル連続性のために以下のスケーリングが必要：

$$
\boxed{\;\frac{f_k^{\rm neq,c}}{f_k^{\rm neq,f}} = \frac{m\,(\tau_c - 1)}{\tau_f - 1}\;}
$$

この係数は Filippova & Hänel (1998) で導出されたもの。物理動粘性連続性と $\delta x_f = \delta x_c/m$ から、$\delta t$ あたりの粘性応力寄与を一致させた結果です。

### 細 → 粗（[lbmblock.c:281-285](../../src/sec4/lbmblock.c#L281-L285)）

粗格子の最下行 $y_c = 0$ は、細格子の $y_f = n_{y,f} - 2$（同じ物理位置）から取得：

$$
f_k^c(2i_c, 0) = f_k^{\rm eq,f}(2i_c, n_{y,f}-2) + \frac{m(\tau_c - 1)}{\tau_f - 1}\!\left[\,f_k^f(2i_c, n_{y,f}-2) - f_k^{\rm eq,f}(2i_c, n_{y,f}-2)\right]
$$

x 方向は細格子の偶数インデックス $i_f = 2 i_c$ が粗格子の $i_c$ に対応。

時間ステップの観点では、粗格子は細格子 $m = 2$ ステップごとに 1 度更新されるため、転送は粗格子の 1 ステップに 1 度だけ実行されます。

### 粗 → 細（[lbmblock.c:380-447](../../src/sec4/lbmblock.c#L380-L447)）

細格子の最上行 $y_f = n_{y,f}$ は、粗格子の $y_c = 1$ から取得。ただし二つの内挿が必要：

#### (a) 時間補間 — Lagrangian 2 次

粗格子は細格子 $m=2$ ステップに 1 回しか進まないため、細格子の **半時刻** $t = n - \tfrac{1}{2}$ における値が必要です。粗格子の過去 3 ステップ $t = n-2, n-1, n$ の値 $\mathrm{ctf}_1, \mathrm{ctf}_2, \mathrm{ctf}_3$ を保存しておき、Lagrange 多項式 $P(t)$ を時刻 $t = n - \tfrac{1}{2}$ で評価：

$$
f_k^{c, n-1/2}(i_c, 1) = \frac{\mathrm{ctf}_3 - 2\,\mathrm{ctf}_2 + \mathrm{ctf}_1}{8} + \frac{\mathrm{ctf}_3 - \mathrm{ctf}_1}{4} + \mathrm{ctf}_2
$$

これは Lagrange 基底を $t = -\tfrac{1}{2}$ で評価したものに等価：

$$
P(-\tfrac{1}{2}) = \tfrac{3}{8}\,\mathrm{ctf}_3 + \tfrac{3}{4}\,\mathrm{ctf}_2 - \tfrac{1}{8}\,\mathrm{ctf}_1
$$

（係数を整理すると上の式と一致）

平衡部 `fc0h` も同様に 3 時刻 `ctf01, ctf02, ctf03` から構築 ([lbmblock.c:391-392](../../src/sec4/lbmblock.c#L391-L392))。$n = 0$（細格子初回ステップ）では 1 ステップ前の値だけ使う簡易版（[lbmblock.c:380-385](../../src/sec4/lbmblock.c#L380-L385)）。

そしてスケーリング則の逆数で非平衡部をスケール：

$$
f_k^f(2i_c, n_{y,f}) = f_k^{\rm eq, c}(i_c, 1; n-1/2) + \frac{\tau_f - 1}{m(\tau_c - 1)}\!\left[\,f_k^{c, n-1/2}(i_c, 1) - f_k^{\rm eq, c}(i_c, 1; n-1/2)\right]
$$

これで **細の偶数 x インデックス** $i_f = 2i_c$ が埋まります。

#### (b) 空間補間 — 自然 3 次スプライン

x 方向の **奇数インデックス** $i_f = 2i_c + 1$ は粗格子に対応点がないため、$n_{x,c}$ 点の値 $a_i = f_k^f(2i_c, n_{y,f})$ から自然 3 次スプライン

$$
S_i(\xi) = a_i + b_i\,\xi + c_i\,\xi^2 + d_i\,\xi^3,\quad \xi \in [0, 1]
$$

を構築し、中点 $\xi = 1/2$ の値を奇数インデックスに代入：

$$
f_k^f(2i_c + 1, n_{y,f}) = a_i + \tfrac{1}{2} b_i + \tfrac{1}{4} c_i + \tfrac{1}{8} d_i
$$

スプライン係数 $c_i$ は連立方程式

$$
\begin{pmatrix}
4 & 1 & & & \\
1 & 4 & 1 & & \\
  & 1 & 4 & 1 & \\
  &   & \ddots & \ddots & \ddots \\
  &   &   & 1 & 4
\end{pmatrix}\!
\begin{pmatrix} c_0 \\ c_1 \\ c_2 \\ \vdots \\ c_{n_{x,c}-1} \end{pmatrix}
=
\begin{pmatrix}
3(a_1 - 2 a_0) \\
3(a_2 - a_1) - 3(a_1 - a_0) \\
3(a_3 - a_2) - 3(a_2 - a_1) \\
\vdots
\end{pmatrix}
$$

を Gauss 消去で解いて求めます（[lbmblock.c:421-436](../../src/sec4/lbmblock.c#L421-L436)）。残りの係数：

$$
b_i = (a_{i+1} - a_i) - \frac{2 c_{i+1} + c_i}{3},\qquad
d_i = \frac{c_{i+1} - c_i}{3}
$$

## 伝播と周期境界

粗・細格子それぞれで [lbmcm.md](lbmcm.md) と同様の伝播を行い、x 方向は周期境界（[lbmblock.c:233, 241, 249, 254, 259, 264](../../src/sec4/lbmblock.c#L233)）：

```c
in = i + 1; if(i == nxc){ in = 0; }   // 右隣
in = i - 1; if(i == 0)  { in = nxc; } // 左隣
```

## マクロ量の評価

[lbmblock.c:288-297, 450-459](../../src/sec4/lbmblock.c#L288-L297) で粗・細それぞれ：

$$
\rho = \sum_k f_k,\qquad
\mathbf{u} = \frac{1}{\rho}\sum_k \mathbf{c}_k\,f_k
$$

## 出力と検証

[lbmblock.c:478-489](../../src/sec4/lbmblock.c#L478-L489) で粗・細を一つの全領域配列 `u[nx][ny]` に統合：

- 下半分 $j \in [0, n_y/2]$ は細格子から $u(i, j) = u_f(2i, 2j)$
- 上半分 $j \in [n_y/2 + 1, n_y]$ は粗格子から $u(i, j) = u_c(i, j - n_y/2 + 1)$

時間 100, 200, 500, 1000, 3000 で中心断面 $u(n_x/2, j)$ をテキストファイル `data{T}` と細格子全断面 `data{T}f` に出力（[lbmblock.c:491-559](../../src/sec4/lbmblock.c#L491-L559)）。

最終的に標準出力に解析解との比較：

```
U[ny    ] : 上壁面（駆動 ut）
U[ny/2+1] : 解析的に ut·(ny/2+1)/ny
U[ny/2  ] : 解析的に (ut+ub)/2 = ut/2
U[ny/2-1] : 解析的に ut·(ny/2-1)/ny
U[0     ] : 下壁面（ub = 0）
```

加えて、中心断面の速度プロファイルが ASCII バーチャートで可視化されます（[lbmblock.c:568-583](../../src/sec4/lbmblock.c#L568-L583)）。

## 計算条件

[lbmblock.c:85, 90, 103-104, 108](../../src/sec4/lbmblock.c#L85-L108)：

| 項目 | 値 |
|---|---|
| 全領域格子 | $n_x \times n_y = 20 \times 32$ |
| 細格子 | $n_{x,f} \times n_{y,f} = 40 \times 32$ |
| 粗格子 | $n_{x,c} \times n_{y,c} = 20 \times 17$ |
| 細分比 | $m = 2$ |
| 粗格子緩和時間 | $\tau_c = 1.2$ |
| 細格子緩和時間 | $\tau_f = 1.9$（自動計算） |
| 動粘性 | $\nu \approx 0.233$（両格子で一致） |
| 上壁速度 | $u_t = 0.01,\ v_t = 0$ |
| 下壁速度 | $u_b = 0,\ v_b = 0$ |
| Reynolds 数 | $Re = u_t\,n_y / \nu \approx 1.37$ |
| 内/外ループ | $30 \times 100 = 3000$ ステップ上限 |
| 収束判定 | $\|\Delta \mathbf{u}\|_\infty < 10^{-10}$ かつ $t > 10000$ |

低 Re（粘性支配）なので、ステップ 100〜3000 で線形 Couette 解に漸近する過渡応答を観察する設計です。

## 実行方法

```powershell
scripts\build_one.cmd src/sec4/lbmblock.c
build\bin\lbmblock.exe
```

出力ファイル（CWD）：

- `data100`, `data200`, `data500`, `data1000`, `data3000` — 中心鉛直線 $u(n_x/2, j)$（全領域、$j = 0\ldots n_y$）
- `data100f`, `data200f`, `data500f`, `data1000f`, `data3000f` — 細格子の中心鉛直線 $u_f(n_{x,f}/2, j_f)$（$j_f = 0\ldots n_{y,f}$）
- `data`, `dataf` — 収束時の最終プロファイル

## 設計判断と注意

- **時間補間が 2 次（Lagrangian quadratic）の理由**: 細格子は粗格子のステップ間に **2 ステップ** 進むため、粗の値が必要なのは「粗の $n-1$ と $n$ の中間 = $n-\tfrac{1}{2}$」の 1 点だけ。2 次補間で十分な精度が出る
- **空間補間が 3 次スプラインの理由**: x 方向は周期境界の連続性が必要。線形補間では C¹ 不連続が応力に影響、3 次スプラインで C² 連続を確保
- **自然境界条件のスプライン**: 始端と終端で $c_0 = ?$ の処理（[lbmblock.c:410-411, 417](../../src/sec4/lbmblock.c#L410-L417)）は本来「2 階微分 = 0」だが、コードは周期境界用ではない簡易実装。Couette の x 一様な流れでは影響は微小
- **スケーリング則の係数 $(\tau - 1)$ vs $(\tau - 1/2)$**: 文献によって $(\tau - 1/2)$ を使う流派もあり、本コードの $(\tau - 1)$ は post-collision 表式に基づく Dupuis-Chopard 系列。テスト結果が解析解に一致する範囲で実用上問題ない
- **オーバーラップ層が 1 セル分のみ**: 粗 $y_c = 0$ と細 $y_f = n_{y,f}-2$ が物理同一点。これより薄いオーバーラップだと衝突→伝播の整合が崩れる
- **DIM = 50** は $\max(n_{x,f}, n_{y,f}) + \alpha = 40 + 10$ の余裕。境界処理 `i+1`, `i-1` のオーバーフロー回避

## 参考

- Filippova, O., Hänel, D. (1998), "Grid Refinement for Lattice-BGK Models", *J. Comput. Phys.*, 147, 219–228 — 多重格子 LBM の標準アルゴリズム原典
- Dupuis, A., Chopard, B. (2003), "Theory and applications of an alternative lattice Boltzmann grid refinement algorithm", *Phys. Rev. E*, 67, 066707 — $(\tau - 1)$ 表式の系列
- Zou, Q., He, X. (1997), "On pressure and velocity boundary conditions for the lattice Boltzmann BGK model", *Phys. Fluids*, 9(6), 1591 — Zou-He 非平衡 bounce-back の原典
- [lbmcm.md](lbmcm.md) — 同じ Couette/cavity 系の単一格子 + 衝突演算子比較
- [cavity.md](cavity.md) — D2Q9 cavity 流れの別実装と Ghia (1982) ベンチマーク
