# sec3: 熱輸送と自然対流の LBM サンプル

本ディレクトリには、双分布関数（double-population）熱格子ボルツマン法に関する 2 本の C プログラムが入っています。両者は D2Q5 thermal LBM の構成を共有しつつ、目的が異なります。

## サンプル一覧

| ファイル | 目的 | ドキュメント |
| --- | --- | --- |
| [lbmnc.c](lbmnc.c) | $Ra = 10^4$ の側面加熱正方キャビティ自然対流（D2Q9 + D2Q5、Boussinesq 浮力、SRT/MRT 切替） | [docs/sec3/lbmnc.md](../../docs/sec3/lbmnc.md) |
| [lbmtherm.c](lbmtherm.c) | 一様流チャネル内 advection-diffusion の D2Q5 ベンチマーク（解析解と比較、サブグリッド壁オフセット $q$） | [docs/sec3/lbmtherm.md](../../docs/sec3/lbmtherm.md) |

## ビルドと実行

リポジトリのルートから次のコマンドで個別実行できます。出力先は既定で `outputs/sec3/<実行ファイル名>/` です。

```powershell
cmd /c scripts\run_one.cmd src\sec3\lbmnc.c
cmd /c scripts\run_one.cmd src\sec3\lbmtherm.c
```

ビルドだけ行いたいときは [scripts/build_one.cmd](../../scripts/build_one.cmd)、まとめてビルドする場合は [scripts/build_all.cmd](../../scripts/build_all.cmd) を使います。

## 解析・可視化スクリプト

| スクリプト | 出力 |
| --- | --- |
| [scripts/plot_lbmnc_schematic.py](../../scripts/plot_lbmnc_schematic.py) | 自然対流解析モデルの模式図 |
| [scripts/plot_lbmnc_results.py](../../scripts/plot_lbmnc_results.py) | 温度・流れ関数・中心線速度・局所 Nusselt 数の 4 + 2 パネル図 |
| [scripts/run_lbmnc_ra_sweep.py](../../scripts/run_lbmnc_ra_sweep.py) | $Ra = 10^3, 10^4, 10^5, 10^6$ の連続実行と de Vahl Davis ベンチマークとの比較 |
| [scripts/plot_lbmtherm_schematic.py](../../scripts/plot_lbmtherm_schematic.py) | lbmtherm.c のチャネル幾何・サブグリッド壁・解析解モードの 3 パネル模式図 |
| [scripts/plot_lbmtherm_results.py](../../scripts/plot_lbmtherm_results.py) | lbmtherm.c の数値解 vs 解析解および誤差場 |
| [scripts/run_lbmtherm_convergence.py](../../scripts/run_lbmtherm_convergence.py) | $n_x \in \{32,48,64,80,96\}$ の格子細分化スイープと $L_2$/$L_\infty$ 収束プロット |

## このディレクトリで扱う物理

- **温度の D2Q5**：5 速度モデル $\mathbf{c}_0 = \mathbf{0},\ \mathbf{c}_{1..4} = (\pm 1, 0), (0, \pm 1)$、重み $w_0 = 1/3,\ w_{1..4} = 1/6$、平衡分布 $g_k^{\mathrm{eq}} = w_k T (1 + 3\mathbf{c}_k\cdot\mathbf{u})$。
- **温度境界条件**：Dirichlet は $g_k^{\mathrm{in}} = -g_{\bar k}^{\mathrm{out}} + 2 w_k T_{\mathrm{wall}}$、Neumann は符号反転なしの bounce-back。サブグリッド壁では $q$ に依存する補間係数で壁面値を再構成（lbmtherm.c。公式は 2 次精度を意図して構築されているが、線形化平衡分布の影響で本コードの実測収束次数は約 1 次）。
- **浮力**：Boussinesq 近似で $F_y = \rho\beta g\,(T - T_{\mathrm{ref}})$、$\rho\beta g$ は $Ra$ から逆算。
- **検証**：lbmnc.c は de Vahl Davis (1983) の正方キャビティベンチマークと比較、lbmtherm.c は閉形式の解析解と直接比較。

両コードを通じて、熱輸送のみの基礎検証（lbmtherm）から、運動量と熱を結合した自然対流の応用（lbmnc）までを一通り体験できる構成になっています。
