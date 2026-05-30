# sec4: 乱流モデルと先端衝突演算子 — cavity, karman, BFS, TG, KH の系統比較

本ディレクトリには、5 つの代表ベンチマーク（lid-driven cavity, backward-facing step, Kármán vortex street, Taylor–Green vortex decay, Kelvin–Helmholtz instability）を、純 BGK と 3 種類の乱流クロージャ（LES Smagorinsky, k-$\varepsilon$, DES）の組合せで解くサンプル群と、先端的な衝突演算子（central-moment LBM）およびマルチブロック格子のサンプルが入っています。

セクション全体は次の 3 軸で構成されています。

1. **問題側**：cavity / BFS / Kármán / Taylor–Green / Kelvin–Helmholtz の 5 種
2. **乱流モデル側**：BGK（無モデル）、LES、k-$\varepsilon$、DES
3. **衝突演算子と格子**：BGK 単一緩和時間、中心モーメント LBM、マルチブロック

## ベース問題（BGK、無モデル）

| ファイル | 物理問題 | ドキュメント |
| --- | --- | --- |
| [cavity.c](cavity.c) | 上壁移動正方キャビティ（Re~384 規模、Ghia ベンチマーク） | [cavity.md](../../docs/sec4/cavity.md) |
| [backward_step.c](backward_step.c) | 後向き段差（reattachment 長さの計測） | [backward_step.md](../../docs/sec4/backward_step.md) |
| [karman.c](karman.c) | 円柱後流の Kármán 渦列（Strouhal 数、スペクトル） | [karman.md](../../docs/sec4/karman.md) |
| [taylor_green.c](taylor_green.c) | 完全周期境界の Taylor–Green 渦減衰 | [taylor_green.md](../../docs/sec4/taylor_green.md) |
| [kelvin_helmholtz.c](kelvin_helmholtz.c) | 速度せん断層の Kelvin–Helmholtz 不安定 | [kelvin_helmholtz.md](../../docs/sec4/kelvin_helmholtz.md) |
| [kelbm.c](kelbm.c) | チャネル流の k-$\varepsilon$ プロトタイプ（壁関数、Dirichlet 境界） | [kelbm.md](../../docs/sec4/kelbm.md) |

## LES バリアント

サブグリッドスケール（SGS）粘性として Smagorinsky モデル $\nu_t = (C_s\Delta)^2 \lvert\overline{S}\rvert$ を加えた D2Q9-BGK 実装です。

| ファイル | ドキュメント |
| --- | --- |
| [cavity_les.c](cavity_les.c) | [cavity_les.md](../../docs/sec4/cavity_les.md) |
| [backward_step_les.c](backward_step_les.c) | [backward_step_les.md](../../docs/sec4/backward_step_les.md) |
| [karman_les.c](karman_les.c) | [karman_les.md](../../docs/sec4/karman_les.md) |
| [taylor_green_les.c](taylor_green_les.c) | [taylor_green_les.md](../../docs/sec4/taylor_green_les.md) |
| [kelvin_helmholtz_les.c](kelvin_helmholtz_les.c) | [kelvin_helmholtz_les.md](../../docs/sec4/kelvin_helmholtz_les.md) |
| [kelbm_les.c](kelbm_les.c) | [kelbm_les.md](../../docs/sec4/kelbm_les.md) |

LES バリアント全体の比較は [les_summary.md](../../docs/sec4/les_summary.md) にまとめてあります。

## k-$\varepsilon$ バリアント

乱流運動エネルギー $k$ と散逸率 $\varepsilon$ をそれぞれ追加の輸送方程式で解く RANS クロージャです。

| ファイル | ドキュメント |
| --- | --- |
| [cavity_keps.c](cavity_keps.c) | （[keps_summary.md](../../docs/sec4/keps_summary.md) で総括） |
| [backward_step_keps.c](backward_step_keps.c) | 同上 |
| [karman_keps.c](karman_keps.c) | 同上 |
| [taylor_green_keps.c](taylor_green_keps.c) | 同上 |
| [kelvin_helmholtz_keps.c](kelvin_helmholtz_keps.c) | 同上 |

## DES バリアント

LES と RANS を領域で切替える Detached Eddy Simulation のハイブリッド実装です。

| ファイル | ドキュメント |
| --- | --- |
| [cavity_des.c](cavity_des.c) | [cavity_des.md](../../docs/sec4/cavity_des.md) |
| [backward_step_des.c](backward_step_des.c) | [backward_step_des.md](../../docs/sec4/backward_step_des.md) |
| [karman_des.c](karman_des.c) | [karman_des.md](../../docs/sec4/karman_des.md) |
| [karman_des_hires.c](karman_des_hires.c) | 高解像度版（同じ DES 実装で格子を密にしたもの） |

ヘッダ [sa_closure.h](sa_closure.h) は Spalart-Allmaras クロージャ係数を共有定義として提供しています。

## 先端衝突演算子・マルチブロック

| ファイル | 目的 | ドキュメント |
| --- | --- | --- |
| [lbmcm.c](lbmcm.c) | 中心モーメント（cascaded）D2Q9 衝突、cavity 流での BGK との精度比較 | [lbmcm.md](../../docs/sec4/lbmcm.md) |
| [lbmcm_karman.c](lbmcm_karman.c) | 中心モーメント LBM による Kármán 渦列、$-5/3$ 慣性領域スペクトルの再現 | [lbmcm_karman.md](../../docs/sec4/lbmcm_karman.md) |
| [lbmblock.c](lbmblock.c) | マルチブロック格子（粗粒/細粒結合）による Couette 流 | [lbmblock.md](../../docs/sec4/lbmblock.md) |

## ビルドと実行

リポジトリのルートから次のコマンドで個別実行できます。出力先は既定で `outputs/sec4/<実行ファイル名>/` です。

```powershell
cmd /c scripts\run_one.cmd src\sec4\karman.c
cmd /c scripts\run_one.cmd src\sec4\lbmcm_karman.c
```

## 解析・可視化スクリプト（抜粋）

| カテゴリ | スクリプト |
| --- | --- |
| Cavity 系 | [plot_cavity_centerline.py](../../scripts/plot_cavity_centerline.py), [plot_cavity_streamlines.py](../../scripts/plot_cavity_streamlines.py) |
| Backward-facing step | [plot_backward_step_streamlines.py](../../scripts/plot_backward_step_streamlines.py), [plot_backward_step_history.py](../../scripts/plot_backward_step_history.py) |
| Kármán 渦列 | [plot_karman_snapshots.py](../../scripts/plot_karman_snapshots.py), [plot_karman_spectrum.py](../../scripts/plot_karman_spectrum.py) |
| Taylor–Green / Kelvin–Helmholtz | [plot_kelvin_helmholtz_snapshots.py](../../scripts/plot_kelvin_helmholtz_snapshots.py), [plot_kelvin_helmholtz_growth.py](../../scripts/plot_kelvin_helmholtz_growth.py) |
| k-$\varepsilon$ チャネル | [plot_kelbm_centerline.py](../../scripts/plot_kelbm_centerline.py), [plot_kelbm_compare_theory.py](../../scripts/plot_kelbm_compare_theory.py), [plot_kelbm_contour.py](../../scripts/plot_kelbm_contour.py), [plot_kelbm_output.py](../../scripts/plot_kelbm_output.py) |
| LES 指標 | [compute_les_metrics.py](../../scripts/compute_les_metrics.py) |
| DES 領域マップ | [plot_des_region_map.py](../../scripts/plot_des_region_map.py) |
| 中心モーメント LBM | [plot_lbmcm_compare.py](../../scripts/plot_lbmcm_compare.py), [plot_lbmcm_distribution.py](../../scripts/plot_lbmcm_distribution.py), [plot_lbmcm_karman_spectrum.py](../../scripts/plot_lbmcm_karman_spectrum.py) |
| マルチブロック | [plot_lbmblock_couette.py](../../scripts/plot_lbmblock_couette.py), [plot_lbmblock_distribution.py](../../scripts/plot_lbmblock_distribution.py) |

## このディレクトリで扱う物理と数値

- **乱流モデル**：BGK で実効粘性を増やすことで乱流をモデル化する SGS 粘性（LES Smagorinsky）、追加の輸送方程式で乱流量を解く RANS（k-$\varepsilon$）、両者を切り替えるハイブリッド DES の 3 方式を、同じベース問題に重ねて比較できます。
- **中心モーメント LBM**：BGK の単一緩和時間に代えて、Galilean-invariant な中心モーメント空間で衝突を行うことで、高 Re での安定性と慣性領域の表現を改善します。
- **マルチブロック格子**：粗粒格子と細粒格子を境界面で接続し、計算量を抑えつつ局所的に高解像度化する手法。
- **代表診断**：Kármán 渦列の Strouhal 数とエネルギースペクトル ($-5/3$ 慣性領域)、Kelvin–Helmholtz 不安定の線形成長率、cavity の中心線速度（Ghia 比較）、backward-facing step の再付着長など、それぞれのベンチマークで一般的な評価指標を計算しています。

サブセクションのテーマ別総括ドキュメントとして [les_summary.md](../../docs/sec4/les_summary.md) と [keps_summary.md](../../docs/sec4/keps_summary.md) が用意されており、5 問題 × 4 モデルの組合せ全体を俯瞰したいときの入口になります。
