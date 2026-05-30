# sec2: 境界条件と基礎流れ — Poiseuille, cavity, 境界条件の比較

本ディレクトリには、LBM の境界条件処理と代表的な 2 次元基礎流れの実装サンプルが入っています。Poiseuille 流（圧力駆動チャネル流）、lid-driven cavity 流、そして 6 種類の境界条件（equilibrium / halfway bounce-back / on-grid / Zou-He / Inamuro / IBL）の系統的比較が含まれ、有限差分法による参考実装（advection 方程式の Upwind / FTCS / Lax-Wendroff / Leap-Frog）も並列に置かれています。

## サンプル一覧

| ファイル | 目的 | ドキュメント |
| --- | --- | --- |
| [fdmadv.c](fdmadv.c) | 1 次元 advection 方程式の有限差分スキーム比較（Upwind, FTCS, Lax-Wendroff, Leap-Frog） | [docs/sec2/fdmadv.md](../../docs/sec2/fdmadv.md) |
| [fdlbm.c](fdlbm.c) | 有限差分版 LBM による Poiseuille 流 | [docs/sec2/fdlbm.md](../../docs/sec2/fdlbm.md) |
| [lbmpoi.c](lbmpoi.c) | 標準 LBM による Poiseuille 流（緩和時間 $\tau$ 依存性、$f^{\mathrm{neq}}$ の評価） | [docs/sec2/lbmpoi.md](../../docs/sec2/lbmpoi.md) |
| [lbmbound.c](lbmbound.c) | 6 種類の境界条件方式を同一テスト問題で系統比較 | [docs/sec2/lbmbound.md](../../docs/sec2/lbmbound.md) |
| [lbmcavi.c](lbmcavi.c) | 上壁移動正方キャビティ（Ghia ベンチマーク、密度変動による compressibility error） | [docs/sec2/lbmcavi.md](../../docs/sec2/lbmcavi.md), [PR サマリ](../../docs/sec2/lbmcavi_pr_summary.md) |

## ビルドと実行

リポジトリのルートから次のコマンドで個別実行できます。出力先は既定で `outputs/sec2/<実行ファイル名>/` です。

```powershell
cmd /c scripts\run_one.cmd src\sec2\lbmcavi.c
cmd /c scripts\run_one.cmd src\sec2\lbmbound.c
```

## 解析・可視化スクリプト

| カテゴリ | スクリプト |
| --- | --- |
| Advection FD 比較 | [plot_fdmadv_courant_compare.py](../../scripts/plot_fdmadv_courant_compare.py) |
| FD-LBM Poiseuille | [plot_fdlbm_profile.py](../../scripts/plot_fdlbm_profile.py) |
| LBM Poiseuille | [plot_lbmpoi_tau_compare.py](../../scripts/plot_lbmpoi_tau_compare.py), [plot_lbmpoi_fneq5_tau_compare.py](../../scripts/plot_lbmpoi_fneq5_tau_compare.py) |
| 境界条件比較 | [plot_lbmbound_all_methods_tau_056.py](../../scripts/plot_lbmbound_all_methods_tau_056.py), [plot_lbmbound_equilibrium_tau_compare.py](../../scripts/plot_lbmbound_equilibrium_tau_compare.py), [plot_lbmbound_halfway_tau_compare.py](../../scripts/plot_lbmbound_halfway_tau_compare.py), [plot_lbmbound_ibl_linear_tau_compare.py](../../scripts/plot_lbmbound_ibl_linear_tau_compare.py), [plot_lbmbound_ibl_quadratic_tau_compare.py](../../scripts/plot_lbmbound_ibl_quadratic_tau_compare.py), [plot_lbmbound_inamuro_tau_compare.py](../../scripts/plot_lbmbound_inamuro_tau_compare.py), [plot_lbmbound_ongrid_tau_compare.py](../../scripts/plot_lbmbound_ongrid_tau_compare.py), [plot_lbmbound_zou_tau_compare.py](../../scripts/plot_lbmbound_zou_tau_compare.py), [plot_lbmbound_nx_study.py](../../scripts/plot_lbmbound_nx_study.py) |
| Lid-driven cavity | [plot_lbmcavi_streamfunction.py](../../scripts/plot_lbmcavi_streamfunction.py), [plot_lbmcavi_ghia_compare.py](../../scripts/plot_lbmcavi_ghia_compare.py), [plot_lbmcavi_ghia_error.py](../../scripts/plot_lbmcavi_ghia_error.py), [plot_lbmcavi_density_change.py](../../scripts/plot_lbmcavi_density_change.py), [compare_lbmcavi_delta_grids.py](../../scripts/compare_lbmcavi_delta_grids.py) |

## このディレクトリで扱う物理と数値

- **緩和時間と粘性**：$\nu = (\tau - 0.5)/3$ の関係を、$\tau$ をスイープしながら Poiseuille 流の放物線速度プロファイルおよび精度で検証。
- **半過程 bounce-back**：壁面が格子間の半セル位置にある前提で、$f_k^{\mathrm{in}}(1,j) = f_{\bar k}^{\mathrm{out}}(0,j)$ の単純な反射を行う最も標準的な no-slip 境界。
- **Zou-He / Inamuro / IBL**：壁面で速度や圧力を Dirichlet 指定したいときに、未知の分布関数を保存則と平衡近似から構成する手法群。IBL は内挿型で、高 $\tau$ 領域での安定性を改善。
- **On-grid bounce-back**：壁面が格子点上にあるとして $f_k$ を反射。半過程との 1 セル分の差が、低次精度・低 $\tau$ で目立つ。
- **Compressibility error**：弱圧縮性 LBM の宿命として、上壁速度 $u_t$ に伴って密度の RMS 変動 $\Delta = \bar\rho^{-1}\sqrt{N^{-1}\sum(\rho - \bar\rho)^2}$ が $O(Ma^2)$ で増大する様子を、cavity 流での測定として確認可能。
- **Ghia ベンチマーク**：Re = 100 の lid-driven cavity 中心線速度を Ghia, Ghia & Shin (1982) の表値と比較し、格子 51×51 程度の粗格子でも全体形状をよく再現することを示します。

このセクションは sec1 の D2Q9-BGK 実装を起点に、境界条件と圧縮性誤差という LBM 特有の論点を一通り体験できる構成です。
