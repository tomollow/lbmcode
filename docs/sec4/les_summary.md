# LES（Smagorinsky）モデルの体系的比較（sec4 6 ケース横断）

## 概要

[src/sec4/](../../src/sec4/) には [keps_summary](keps_summary.md) で扱った 6 ケースに対し、標準型 Smagorinsky LES を組み合わせた `*_les.c` を追加実装しました。本ドキュメントは**同じ Smagorinsky SGS モデルを異なる流れに適用したときの挙動**を一望し、k-ε との対照を通じて「Smagorinsky は局所応答型 SGS として何を modeling するか」を定量的に整理するクロスケース比較です。

SGS モデルは標準 Smagorinsky を採用：

$$ \nu_t = (C_s\,\Delta)^2 \sqrt{2\,S_{ij}S_{ij}}, \quad C_s = 0.16,\ \Delta = 1 \text{ LU} $$

van Driest 減衰や動的係数決定は導入していません。$\nu_t$ は BGK 衝突の局所緩和時間 $\tau_{\rm eff} = 1/2 + 3(\nu_0 + \nu_t)$ に取り込まれます。

## 比較表（LES の効果量）

| ケース | 物理レジーム | $Re$ 系 | $\nu_t/\nu_0$ (平均) | LES による主要量への影響 | 物理的妥当性 |
|---|---|---:|---:|---|---|
| [kelbm](kelbm.md) | 低 Re チャンネル流れ | $Re_\tau$ ≈ 22 | 0.005 | $u_{\max}$ 層流予測比 **0.74**（k-ε と同等：いずれも未収束過渡） | ❓ 過渡同等 |
| [taylor_green](taylor_green.md) | 等方減衰渦（周期、外力なし）| $k=2\pi/L$ | 9e-5 | 解析減衰率の **1.00 倍**（解析解を保存）| ✅ 妥当 |
| [kelvin_helmholtz](kelvin_helmholtz.md) | 遷移せん断不安定 | $Re_{\rm shear}$ ≈ 60 | 0.006 | ピーク振幅 **0.93 倍**（7% 抑制）| ✅ 妥当 |
| [cavity](cavity.md) | 層流蓋駆動再循環 | $Re$ ≈ 384 | 0.001 | $\psi_{\min}$ **1.00 倍**（実質無作用）| ✅ 妥当 |
| [backward_step](backward_step.md) | 層流ステップ後流 | $Re_H$ ≈ 56 | 0.001 | 再付着長 $x_R/H$ **1.00 倍**（変化なし）| ✅ 妥当 |
| [karman](karman.md) | 周期渦放出（Hopf 分岐後）| $Re_D$ ≈ 126 | 0.004 | 振幅 **0.86 倍**（14% 抑制）、周波数は同じ | ✅ 妥当 |

すべての $\nu_t/\nu_0$ は LES 履歴の時間平均（バルク平均）です。最大値も $0.01$ 未満で、分子粘性に対する $\nu_t$ の寄与はどのケースでも $1\%$ 程度に留まります。

## 観察される傾向

### 1. Smagorinsky はほぼ「眠った」状態に留まる

全 6 ケースで $\nu_t/\nu_0 \lesssim 0.01$（k-ε の $0.04$–$0.16$ より 1–2 桁小さい）。これは:

$$ \nu_t = (C_s\Delta)^2\,|S| \approx 0.026 \cdot |S|_{\rm LU} $$

において LBM の格子スケール $\Delta = 1$ で $|S|$ が $|U|/L_{\rm characteristic} \lesssim 10^{-2}$ 程度になる低 Re レジームでは $\nu_t$ が極小に留まるためです。**SGS としての本来の役割**——格子で解像できない散逸の補完——は、 Re が大幅に高くなり LBM が under-resolved となるまで活性化しません。

### 2. 遷移流れ・不安定でも LES は控えめに振る舞う

KH (7%)、Karman (14%) では確かに振幅が減衰しますが、k-ε の対応値（21%、86%）と比べて 1 桁小さい：

- **KH**: ピーク振幅は k-ε `0.79 倍` に対し LES `0.93 倍` — 渦の roll-up はほぼ DNS と区別がつかない
- **Karman**: 振幅は k-ε `0.14 倍`（86% 抑制）に対し LES `0.86 倍`（14% 抑制） — 渦放出が物理的に残る
- **理由**: Smagorinsky は瞬時の $|S|$ に応答するため、渦輸送による $k$ の蓄積が起きない。**輸送型 RANS （k-ε）が時間積分した「過去の」乱流エネルギーで現在の流れを抑制する**のに対し、**Smagorinsky は瞬時の歪み速度に局所的に応答する**

### 3. 等方減衰（Taylor-Green）では解析解を保存

Taylor-Green では k-ε が解析減衰の **0.75 倍**（25% 加速減衰）に対し、LES は **1.00 倍**（解析解と区別がつかない）。これは:

- LES 平均 $\nu_t/\nu_0 \approx 9 \times 10^{-5}$ — 実質ゼロ
- 渦が減衰するとともに $|S|$ も小さくなり、$\nu_t$ がさらに減る自己終息的挙動
- 一方 k-ε では一度生成された $k$ が散逸 $\varepsilon$ まで時間遅れで残り、能動的に減衰を加速

この差は標準 k-ε が形式的にせん断流向け RANS であり等方乱流に不適である**理論的限界が定量化された結果**と言えます。

### 4. 壁関数なしの代償

Smagorinsky には壁関数機構がない（van Driest 減衰を入れていないため）。これは:

- **長所**: 設定パラメータが $C_s$ 1 つで済む、複雑な壁処理不要
- **短所**: 壁第1セルで $\nu_t$ がゼロにならない — 壁せん断力を **過小評価**しがち
- 本実装の低 Re ケース（kelbm/cavity/step/karman）では $|S|$ が壁近傍でも小さいため実害は軽微

## LES と k-ε の効果スケーリング比較

| ケース | $\nu_t/\nu_0$ (LES) | $\nu_t/\nu_0$ (k-ε) | LES 抑制 | k-ε 抑制 | k-ε / LES 比 |
|---|---:|---:|---:|---:|---:|
| Cavity | 0.001 | 0.04 | ~0% | 1% | — |
| Step | 0.001 | 0.034 | ~0% | 2% | — |
| TG | 9e-5 | 0.11 | ~0% | 25% | ∞ |
| KH | 0.006 | 0.05 | 7% | 21% | 3.0× |
| Karman | 0.004 | 0.048 | 14% | 86% | 6.1× |
| kelbm | 0.005 | 0.16 | (26%, 過渡) | (26%, 過渡) | 1.0× |

**k-ε / LES 比** は同じ流れに対する k-ε の追加散逸の大きさを示します。遷移/周期不安定（KH, Karman）で 3–6 倍、等方減衰では無限大に発散——これは k-ε が**輸送方程式で過去の歪みの履歴を保持する**ため、瞬時応答型の Smagorinsky と比較して持続的な散逸源となるためです。

kelbm のみ k-ε / LES 比が 1.0 となる理由：[kelbm.c](../../src/sec4/kelbm.c) は k-ε 輸送方程式を解いていますが、BGK 衝突の $\tau$ には $\nu_t$ が戻されておらず（kelbm.c の `stream_collide` は定数 OMEGA を使用）、`u_max/u_lam ≈ 0.74` は単に NSTEPS=30000 が層流定常時間 $T \sim H^2/\nu_0 \approx 2.2 \times 10^5$ ステップに対し未収束であるための過渡値です。LES 版（[kelbm_les.c](../../src/sec4/kelbm_les.c) は $\tau_{\rm eff}$ で結合）でも同じ過渡値が出ているのは、$\nu_t$ がきわめて小さく結合の有無が結果に効かないためです。

## 設計上の共通パターン

k-ε 版との差分のみ列挙：

### コア構造
- **状態変数の削減**: $k$, $\varepsilon$ 配列なし。$\nu_t$ 場 `nut_field[]` のみ
- **輸送方程式なし**: `update_les()` は単に勾配計算と $\nu_t = C_s^2 |S|$ を実行する 1 パス
- **時間ステップ制約緩和**: k-ε で必要だった `KEPS_DT = 0.05` の安定性制限が消える — LES は陽的代数式のみ

### BGK 結合
- **同一の局所 $\tau_{\rm eff}$**: `tau_eff = 0.5 + 3.0 * (nu0 + nut_field[i])` は k-ε 版と完全同一
- **壁関数なし**: `apply_wall_function` を持たない（Smagorinsky 単体には壁処理機構がない）

### CSV 出力
- スナップショット列から `k, eps` を除外、`nut` のみを残す
- 履歴ファイルからも `k_mean, eps_mean` を除外
- ファイル名は `*_keps_*` を `*_les_*` に置換

### 実装行数の縮小

| ファイル | keps 行数 | les 行数 | 削減率 |
|---|---:|---:|---:|
| kelbm | 243 (kelbm.c) | 156 | 36% |
| taylor_green | 242 | 199 | 18% |
| kelvin_helmholtz | 259 | 209 | 19% |
| cavity | 315 | 222 | 30% |
| backward_step | 346 | 244 | 29% |
| karman | 322 | 234 | 27% |

平均 27% 減 — 輸送方程式、壁関数、初期 seed 計算、`k_new/eps_new` バッファが不要になるため。

## 物理的限界

Smagorinsky 単体の既知の限界：

- **過剰散逸近傍壁**: 壁第1セルで $\nu_t$ が消えない。本実装の低 Re では $|S|$ 自体が小さく問題顕在化せず
- **遷移現象に鈍感**: 速度勾配のみに応答するため、層流から乱流への遷移メカニズム（不安定モードの選択的増幅）には介入しない
- **$C_s$ は流れに依存**: 標準 0.16 は等方乱流向け。せん断流（KH, kelbm）では 0.1 程度が推奨
- **格子に依存**: $\Delta = $ 格子間隔と仮定。AMR や非均一格子では再定義必要

これらを補正するには：

- **van Driest 減衰**: $C_s \to C_s [1 - \exp(-y^+/A^+)]$ で壁近傍を抑制（$A^+ \approx 26$）
- **WALE モデル**: 壁近傍で自然に $\nu_t \to 0$
- **動的 Smagorinsky**: $C_s$ をフィルタ比から決定（Germano et al. 1991）

## 教育的ポジショニング

LES シリーズの学習ゴール（k-ε シリーズと相補的）：

1. **SGS モデルの局所性**: 輸送なしの瞬時応答という設計思想（taylor_green, kelvin_helmholtz の対照）
2. **低 Re での「眠った」挙動**: $|S|$ がそもそも小さい流れではモデルが効かない（cavity, step）
3. **過剰散逸を避ける重要性**: 渦放出（karman）で SGS が時間積分しない利点
4. **モデル選択の指針**: 同じ流れに RANS と LES を当てて比較できる稀有な教材構成（[keps_summary](keps_summary.md) と本 doc）

「LES = 高 Re で初めて意味がある」「k-ε = 平均流に強いが瞬時構造を smear する」という**理論的特性が、6 ケースのコード実装で具体的な数値として再現される**ことを体験できます。

## 検証メトリクスの再現

本 doc の表値を再生成するには：

```powershell
# 各ランナーは -LesOnly で LES 版だけビルド・実行・プロットを行う
.\scripts\run_kelbm.ps1            -LesOnly
.\scripts\run_taylor_green.ps1     -LesOnly
.\scripts\run_kelvin_helmholtz.ps1 -LesOnly
.\scripts\run_cavity.ps1           -LesOnly
.\scripts\run_backward_step.ps1    -LesOnly
.\scripts\run_karman.ps1           -LesOnly

# pure / k-eps / LES の比較表値を算出
.\.venv\Scripts\python.exe scripts\compute_les_metrics.py
```

ランナーは `-PureOnly` / `-KepsOnly` / `-LesOnly` のいずれか 1 つを受け付け、無指定なら全 variant を実行します。

[scripts/compute_les_metrics.py](../../scripts/compute_les_metrics.py) は LES / pure / k-ε の結果を読んで上の比較表の数値を算出します。

## 参考文献

- Smagorinsky 原典: Smagorinsky (1963) "General circulation experiments with the primitive equations", Monthly Weather Review
- LES 入門: Sagaut *Large Eddy Simulation for Incompressible Flows* (2006)
- LBM × LES: Yu, Girimaji & Luo (2005) "DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method", J. Comp. Phys.
- Dynamic Smagorinsky: Germano, Piomelli, Moin & Cabot (1991) "A dynamic subgrid-scale eddy viscosity model", Phys. Fluids A
- WALE モデル: Nicoud & Ducros (1999) "Subgrid-scale stress modelling based on the square of the velocity gradient tensor", Flow Turbul. Combust.
