# 図生成スクリプト一覧

このリポジトリでは、説明用の図や比較プロットを `scripts` 配下の Python スクリプトで再生成できます。ビルドや実行を伴うものは、各スクリプトが必要な処理を内部で行います。

## 実行方法

リポジトリのルートで、次のように実行します。

```powershell
d:/work/LBMcode/.venv/Scripts/python.exe scripts/<script-name>.py
```

## 現在ある図生成スクリプト

### sec1

- [scripts/plot_lbmtv_nx_study.py](../scripts/plot_lbmtv_nx_study.py)
  - 用途: Taylor vortex の格子点数 `nx` と誤差 `Erru` の関係を集計して図にする
  - 主な出力: [docs/assets/sec1/lbmtv_nx_errors.png](assets/sec1/lbmtv_nx_errors.png)
  - 補助出力: [docs/sec1/generated/lbmtv_nx_errors.md](sec1/generated/lbmtv_nx_errors.md), [docs/sec1/generated/lbmtv_nx_errors.csv](sec1/generated/lbmtv_nx_errors.csv)

### sec2

- [scripts/plot_poiseuille_comparison.py](../scripts/plot_poiseuille_comparison.py)
  - 用途: Poiseuille flow の比較図を作る
  - 主な出力: [docs/assets/sec2/poiseuille_profile_comparison.png](assets/sec2/poiseuille_profile_comparison.png)

- [scripts/plot_fdlbm_profile.py](../scripts/plot_fdlbm_profile.py)
  - 用途: `fdlbm.c` の中心断面速度分布と解析解の比較図を作る
  - 主な出力: [docs/assets/sec2/fdlbm_profile.png](assets/sec2/fdlbm_profile.png)

- [scripts/plot_lbmpoi_tau_compare.py](../scripts/plot_lbmpoi_tau_compare.py)
  - 用途: `tau` を変えたときの Poiseuille flow 速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmpoi_tau_compare.png](assets/sec2/lbmpoi_tau_compare.png)

- [scripts/plot_lbmbound_all_methods_tau_056.py](../scripts/plot_lbmbound_all_methods_tau_056.py)
  - 用途: `tau = 0.56` で 7 つの境界条件の速度分布を 1 枚にまとめた比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_all_methods_tau_056.png](assets/sec2/lbmbound_all_methods_tau_056.png)

- [scripts/plot_lbmbound_nx_study.py](../scripts/plot_lbmbound_nx_study.py)
  - 用途: `tau = 0.56` で 7 つの境界条件の格子収束を集計し、表と図にまとめる
  - 主な出力: [docs/assets/sec2/lbmbound_nx_errors.png](assets/sec2/lbmbound_nx_errors.png)
  - 補助出力: [docs/sec2/generated/lbmbound_nx_errors.md](sec2/generated/lbmbound_nx_errors.md), [docs/sec2/generated/lbmbound_nx_errors.csv](sec2/generated/lbmbound_nx_errors.csv)

- [scripts/plot_lbmbound_equilibrium_tau_compare.py](../scripts/plot_lbmbound_equilibrium_tau_compare.py)
  - 用途: Equilibrium 境界条件で `tau` を変えたときの速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_equilibrium_tau_compare.png](assets/sec2/lbmbound_equilibrium_tau_compare.png)

- [scripts/plot_lbmbound_ongrid_tau_compare.py](../scripts/plot_lbmbound_ongrid_tau_compare.py)
  - 用途: On-grid bounce back 境界で `tau` を変えたときの速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_ongrid_tau_compare.png](assets/sec2/lbmbound_ongrid_tau_compare.png)

- [scripts/plot_lbmbound_inamuro_tau_compare.py](../scripts/plot_lbmbound_inamuro_tau_compare.py)
  - 用途: Inamuro の滑りなし境界条件で `tau` を変えたときの速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_inamuro_tau_compare.png](assets/sec2/lbmbound_inamuro_tau_compare.png)

- [scripts/plot_lbmbound_halfway_tau_compare.py](../scripts/plot_lbmbound_halfway_tau_compare.py)
  - 用途: Half-way bounce back 境界で `tau` を変えたときの速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_halfway_tau_compare.png](assets/sec2/lbmbound_halfway_tau_compare.png)

- [scripts/plot_lbmbound_ibl_linear_tau_compare.py](../scripts/plot_lbmbound_ibl_linear_tau_compare.py)
  - 用途: 線形補間バウンスバック境界で `tau` を変えたときの速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_ibl_linear_tau_compare.png](assets/sec2/lbmbound_ibl_linear_tau_compare.png)

- [scripts/plot_lbmbound_ibl_quadratic_tau_compare.py](../scripts/plot_lbmbound_ibl_quadratic_tau_compare.py)
  - 用途: 二次補間バウンスバック境界で `tau` を変えたときの速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_ibl_quadratic_tau_compare.png](assets/sec2/lbmbound_ibl_quadratic_tau_compare.png)

- [scripts/plot_lbmbound_zou_tau_compare.py](../scripts/plot_lbmbound_zou_tau_compare.py)
  - 用途: Zou-He の非平衡バウンスバック境界で `tau` を変えたときの速度分布比較図を作る
  - 主な出力: [docs/assets/sec2/lbmbound_zou_tau_compare.png](assets/sec2/lbmbound_zou_tau_compare.png)

- [scripts/plot_non_equilibrium_distribution.py](../scripts/plot_non_equilibrium_distribution.py)
  - 用途: 非平衡 bounce-back の模式図を作る
  - 主な出力: [docs/assets/sec2/non_equilibrium_distribution.png](assets/sec2/non_equilibrium_distribution.png)

- [scripts/plot_lbmpoi_fneq5_tau_compare.py](../scripts/plot_lbmpoi_fneq5_tau_compare.py)
  - 用途: 非平衡分布関数 $f_5^{neq}$ の `tau` 比較図を作る
  - 主な出力: [docs/assets/sec2/lbmpoi_fneq5_tau_compare.png](assets/sec2/lbmpoi_fneq5_tau_compare.png)

- [scripts/plot_lbmcavi_streamfunction.py](../scripts/plot_lbmcavi_streamfunction.py)
  - 用途: `lbmcavi.c` の流れ関数等値線図と中心線速度プロファイルを作る
  - 前提: [outputs/sec2/lbmcavi](../outputs/sec2/lbmcavi) に `datau`, `datav`, `datas` があること
  - 主な出力: [docs/assets/sec2/lbmcavi_streamfunction.png](assets/sec2/lbmcavi_streamfunction.png)
  - 補助出力: [docs/sec2/generated/lbmcavi_ghia_re100_comparison.csv](sec2/generated/lbmcavi_ghia_re100_comparison.csv)

- [scripts/plot_lbmcavi_ghia_compare.py](../scripts/plot_lbmcavi_ghia_compare.py)
  - 用途: `lbmcavi.c` の中心線速度を Ghia (1982) の Re = 100 ベンチマークと比較する
  - 前提: [outputs/sec2/lbmcavi](../outputs/sec2/lbmcavi) に `datau`, `datav` があること
  - 主な出力: [docs/assets/sec2/lbmcavi_ghia_compare.png](assets/sec2/lbmcavi_ghia_compare.png)

- [scripts/plot_lbmcavi_ghia_error.py](../scripts/plot_lbmcavi_ghia_error.py)
  - 用途: `lbmcavi.c` の Ghia 比較表から中心線速度の絶対誤差を可視化する
  - 前提: [docs/sec2/generated/lbmcavi_ghia_re100_comparison.csv](sec2/generated/lbmcavi_ghia_re100_comparison.csv) があること
  - 主な出力: [docs/assets/sec2/lbmcavi_ghia_error.png](assets/sec2/lbmcavi_ghia_error.png)

- [scripts/plot_lbmcavi_density_change.py](../scripts/plot_lbmcavi_density_change.py)
  - 用途: 上壁速度 `u_x^w` に対する平均密度変化 `Δ = (1/ρ̄)√((1/N)Σ(ρ-ρ̄)^2)` の関係を作図する
  - 前提: `lbmcavi.c` が上壁速度をコマンドライン引数で受け取れること
  - 主な出力: [docs/assets/sec2/lbmcavi_density_change_vs_uwx.png](assets/sec2/lbmcavi_density_change_vs_uwx.png)
  - 補助出力: [docs/sec2/generated/lbmcavi_density_change_vs_uwx.csv](sec2/generated/lbmcavi_density_change_vs_uwx.csv)

## 補足

- スクリプトによっては C 実行ファイルのビルドを先に行います。
- 実行時の生データは主に `outputs/` 配下に残し、コミット対象にしたい表や図は `docs/` 配下へ複製しています。
- `lbmtv_nx_errors.png` の具体的な生成手順は [docs/sec1/lbmtv.md](sec1/lbmtv.md) にも記載しています。