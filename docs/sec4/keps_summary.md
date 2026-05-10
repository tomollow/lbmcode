# k-ε モデルの体系的比較（sec4 6 ケース横断）

## 概要

[src/sec4/](../../src/sec4/) には LBM + 標準 $k$-$\varepsilon$ 乱流モデルの実装が 6 ケース揃っており、それぞれが異なる物理レジーム（層流・遷移・不安定・周期渦放出など）を扱います。本ドキュメントは**同じ k-ε モデルを異なる流れに適用したときの挙動**を一望し、「k-ε はどんな流れに有効で、どんな流れに過剰散逸となるか」を定量的に整理するクロスケース比較です。

## 比較表（k-ε の効果量）

| ケース | 物理レジーム | $Re$ 系 | $\nu_t/\nu_0$ | k-ε による主要量への影響 | 物理的妥当性 |
|---|---|---:|---:|---|---|
| [kelbm](kelbm.md) | 低 Re チャンネル流れ（壁関数あり）| $Re_\tau$ ≈ 22 | 0.16 | $u_{\max}$ を層流予測の **0.74 倍** に抑制 | ✅ 良性 |
| [taylor_green](taylor_green.md) | 等方減衰渦（周期、外力なし）| $k=2\pi/L$ | 0.11 | 解析減衰率の **0.75 倍**（25% 加速減衰）| ⚠️ RANS 範囲外 |
| [kelvin_helmholtz](kelvin_helmholtz.md) | 遷移せん断不安定 | $Re_{\rm shear}$ ≈ 60 | 0.05 | ピーク振幅 **0.79 倍**（21% 抑制）| ⚠️ 過剰散逸 |
| [cavity](cavity.md) | 層流蓋駆動再循環 | $Re$ ≈ 384 | 0.04 | $\psi_{\min}$ **1.01 倍**（実質無作用）| ✅ 良性 |
| [backward_step](backward_step.md) | 層流ステップ後流 | $Re_H$ ≈ 56 | 0.034 | 再付着長 $x_R/H$ **0.97 倍**（2% 短縮）| ✅ 良性 |
| [karman](karman.md) | 周期渦放出（Hopf 分岐後）| $Re_D$ ≈ 127 | 0.048 | 振幅 **0.144 倍**（86% 抑制）、周波数は同じ | ⚠️ 過剰散逸（最も顕著）|

## 観察される傾向

### 1. 層流レジームでは k-ε はほぼ無作用

Cavity・Step・低 Re kelbm のような**定常または準定常な層流**では、ひずみ速度 $S^2$ が小さく $P_k = \nu_t S^2$ が dissipation $\varepsilon$ と釣り合った位置で $\nu_t$ が小さく安定。kelbm では壁関数で $k$ を強制注入しているため $\nu_t/\nu_0$ がやや大きい（0.16）が、それでも流れを破綻させずに自然に整合。

### 2. 遷移流れ・不安定では過剰散逸

KH の渦巻き形成、Karman の渦放出のような**遷移または周期不安定**では、k-ε の渦粘性が瞬時の渦構造を smearing する方向に働く：

- **KH**: 渦の roll-up はほぼ同じ（モード 4 の 4 渦）だが、ピーク振幅が下がり、後期の渦合体・散逸が早まる
- **Karman**: 渦放出周波数（Strouhal $\sim 0.14$）は変わらない（時間平均構造として）が、振幅が **86%** 抑制 — 最も劇的

これは標準 k-ε が**時間平均流向け RANS 閉鎖**であり、瞬時の渦構造の捕捉に向かないという理論通りの挙動。

### 3. 等方減衰（Taylor-Green）は中間

完全周期・外力なし・解析解が存在する Taylor-Green では、k-ε による $\nu_t$ が分子粘性に上乗せされ、減衰が**指数的に加速**（25%）。形式的には標準 k-ε はせん断流向けで TG には不適だが、ソース項の挙動チェックとしては機能。

## k-ε 効果のスケーリング

$\nu_t/\nu_0$ の値だけでなく、**「速度抑制への変換効率」**がレジームで異なる：

| ケース | $\nu_t/\nu_0$ | 主要量抑制 | 効率（抑制 / $\nu_t/\nu_0$） |
|---|---:|---:|---:|
| Cavity | 0.04 | 1% | 0.25 |
| Step | 0.034 | 2% | 0.6 |
| TG | 0.11 | 25% | 2.3 |
| KH | 0.05 | 21% | 4.2 |
| **Karman** | **0.048** | **86%** | **17.9** |

Karman は同じ $\nu_t/\nu_0$ レベルでも極端に大きな抑制を見せる — これは渦放出が**振動現象**で、平均流ではなく振動振幅を減衰させる k-ε の特性が劇的に現れるため。

## 設計上の共通パターン

5 ケース（+ kelbm）で再利用された実装パターン：

### LBM コア
- **D2Q9 BGK 衝突**: 全ケース共通の `feq = w * rho * (1 + 3 c·u + 4.5 (c·u)² - 1.5 |u|²)`
- **f / f2 ポインタスワップ**: 各ステップで `double *tmp = f; f = f2; f2 = tmp;` で 2× 高速化
- **`%.9g` 精度の CSV 出力**: 再現性のため

### 境界条件
- **halfway bounce-back**（kelbm, KH, cavity, step, karman）: 反射時 `f2[i*NDIR + opp[d]] = post`
- **Ladd 動壁補正**（cavity の上端）: `post - 6 w[d] rho c_x U_lid`
- **`solid[]` マスク**（step, karman）: 任意形状の障害物

### 駆動と初期条件
- **Guo forcing**: チャンネル流れ（kelbm, step, karman）に体積力 $F_x$、$\tau_{\rm eff}$ 付き $F_k$ 補正
- **Asymmetric perturbation**: KH のモード 4 sin、Karman の wake Gaussian — 不安定モード励起

### k-ε 結合
- **局所 $\tau_{\rm eff} = 1/2 + 3(\nu_0 + \nu_t)$**: BGK 衝突に渦粘性を取り込む
- **Dirichlet 壁関数**: $k_{\rm wall} = u_\tau^2/\sqrt{C_\mu}$、$\varepsilon_{\rm wall} = u_\tau^3/(\kappa\,\Delta y)$
- **`apply_wall_function` の max ロジック**: コーナーセルが両壁の寄与を取れるよう reset しない
- **対流項は 1 次風上、拡散項は 2 次中心**: 安定性とコストのバランス

### Runner / Plot
- **`scripts/run_<name>.ps1`**: build → run pure / k-eps → plot 生成、`-PureOnly`/`-KepsOnly`/`-SkipPlot` 対応、両者同時指定の guard あり
- **Path ベース Python スクリプト**: `ROOT_DIR = Path(__file__).resolve().parents[1]`、`outputs/sec4/` と `docs/assets/sec4/` の両方に保存
- **ベクトル化 pivot**: `df.pivot(index='y', columns='x', values=...).reindex(...)` で `iterrows` を回避

## 物理的限界（共通）

すべての k-ε 版で：

- **本来の RANS は時間平均流向け**で、瞬時の渦放出（Karman）や周期不安定（KH）の精密な捕捉には不適
- **壁関数は形式的に $y^+ \gtrsim 30$ で有効**だが、本実装の典型 $y^+_1 \approx 0.5$–2 は対数領域外。実態は「$k$ 注入器」として機能
- **BGK は $\tau \to 0.5$ で不安定化** — 高 Re 域には MRT/regularized 衝突演算子推奨
- **2D 限定**: Karman は $Re > 200$ で 3D 不安定（Mode A/B）が現れるため上限。完全乱流チャンネルや実用 BFS は 3D LES が必要

## 教育的ポジショニング

このシリーズは **LBM + RANS の入門教材**として設計されており、以下の学習ゴールに対応：

1. **LBM の基本実装**: BGK 衝突、平衡分布、halfway bounce-back の動作確認（kelbm, TG）
2. **境界条件の多様性**: 周期、固定壁、動壁、固体ブロック、円柱（曲面 staircase）（cavity → karman）
3. **流体不安定の再現**: 線形理論との対比（KH, Karman の Strouhal）
4. **乱流モデルの限界の理解**: k-ε の良性/過剰散逸レジームの定量比較（本まとめ doc）

研究レベルの精度には不十分だが、**「k-ε モデルが何を modeling しているか」** を実体験で理解する素材としては唯一無二の構成です。

## 参考文献まとめ

各ケース doc に詳細：

- LBM 標準: He & Luo (1997), Krüger et al. *The Lattice Boltzmann Method* (2017)
- k-ε モデル: Launder & Spalding (1974), Pope *Turbulent Flows* (2000)
- Taylor-Green: Taylor & Green (1937), Brachet et al. (1983)
- Kelvin-Helmholtz: Drazin & Reid *Hydrodynamic Stability* (1981)
- Cavity: Ghia, Ghia & Shin (1982)
- BFS: Armaly et al. (1983)
- Kármán: Williamson (1996), Roshko (1954)
