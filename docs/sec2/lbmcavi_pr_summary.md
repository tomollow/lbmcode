# lbmcavi 更新記録

## 概要

`lbmcavi.c` まわりの実行系、説明ドキュメント、図版生成スクリプト、比較用 CSV を追加・整理した。

## 変更内容

- `lbmcavi.c` がコマンドライン引数で `u_t` と `nx` を受け取れるように変更
- `lbmcavi.c` の配列を静的領域へ移し、`nx = 256` を扱えるように `DIM` を拡張
- `error` 出力を `(u_t, rho_avg, delta_rms, err)` の 4 列へ整理
- `run_one.cmd` の出力先を `outputs/sec2/lbmcavi` へ安定化
- `docs/sec2/lbmcavi.md` を新規追加し、数式・境界条件・実行例・出力説明を整理
- 流れ関数図、Ghia 比較図、Ghia 誤差図、密度変化図の生成スクリプトを追加
- 公開用画像を `docs/assets/sec2/` に追加
- Ghia 比較 CSV、密度変化 CSV、格子比較 CSV を `docs/sec2/generated/` に追加
- 参考 PDF `docs/assets/sec2/9401003v1.pdf` を追加
- `README.md` と `docs/plot_generation.md` に `lbmcavi` 関連項目を追記

## 主な追加ファイル

- [docs/sec2/lbmcavi.md](lbmcavi.md)
- [scripts/plot_lbmcavi_streamfunction.py](../../scripts/plot_lbmcavi_streamfunction.py)
- [scripts/plot_lbmcavi_ghia_compare.py](../../scripts/plot_lbmcavi_ghia_compare.py)
- [scripts/plot_lbmcavi_ghia_error.py](../../scripts/plot_lbmcavi_ghia_error.py)
- [scripts/plot_lbmcavi_density_change.py](../../scripts/plot_lbmcavi_density_change.py)
- [scripts/compare_lbmcavi_delta_grids.py](../../scripts/compare_lbmcavi_delta_grids.py)

## 生成物

- [docs/assets/sec2/lbmcavi_streamfunction.png](../assets/sec2/lbmcavi_streamfunction.png)
- [docs/assets/sec2/lbmcavi_ghia_compare.png](../assets/sec2/lbmcavi_ghia_compare.png)
- [docs/assets/sec2/lbmcavi_ghia_error.png](../assets/sec2/lbmcavi_ghia_error.png)
- [docs/assets/sec2/lbmcavi_density_change_vs_uwx.png](../assets/sec2/lbmcavi_density_change_vs_uwx.png)
- [docs/sec2/generated/lbmcavi_ghia_re100_comparison.csv](generated/lbmcavi_ghia_re100_comparison.csv)
- [docs/sec2/generated/lbmcavi_density_change_vs_uwx.csv](generated/lbmcavi_density_change_vs_uwx.csv)
- [docs/sec2/generated/lbmcavi_delta_grid_comparison.csv](generated/lbmcavi_delta_grid_comparison.csv)

## 関連コミット

- `d77c618` Make lbmcavi configurable and fix run output path
- `88cedab` Add lbmcavi documentation, plots, and analysis scripts
- `3862c83` Add cavity reference PDF asset

## 確認内容

- `docs/sec2/lbmcavi.md` に Markdown エラーがないことを確認
- `scripts/plot_lbmcavi_density_change.py` の実行により `nx = 256` の密度変化図と CSV を再生成
- `origin/main` へ push 済み