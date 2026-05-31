# sec6: 没入境界 LBM — 任意形状物体まわりの流れ

本ディレクトリには、Immersed Boundary–Lattice Boltzmann Method (IB-LBM) のサンプルが入っています。流体格子とは独立に Lagrangian 点で表現された物体境界に、適切な体積力を加えることで no-slip 条件を満たす方法群です。3 種類の衝突演算子（SRT, MRT, TRT）と 2 種類の体積力法（Direct Forcing, Implicit Correction）の組合せ、および落下粒子のサンプルが揃っています。

## サンプル一覧

### 円筒 Couette 流（解析解との比較）

内外 2 重円筒の Couette 流（外側静止、内側回転）を IB-LBM で解き、回転 Stokes 流の解析解と比較するベンチマークです。

| ファイル | 体積力法 | 衝突演算子 | 解説 |
| --- | --- | --- | --- |
| [iblbm2cdfSRT.c](iblbm2cdfSRT.c) | Direct Forcing | SRT | [docs/sec6/iblbm2cdfSRT.md](../../docs/sec6/iblbm2cdfSRT.md) |
| [iblbm2cicMRT.c](iblbm2cicMRT.c) | Implicit Correction | MRT | [docs/sec6/iblbm2cicMRT.md](../../docs/sec6/iblbm2cicMRT.md) |
| [iblbm2cicTRT.c](iblbm2cicTRT.c) | Implicit Correction | TRT | [docs/sec6/iblbm2cicTRT.md](../../docs/sec6/iblbm2cicTRT.md) |

円筒 Couette 流ベンチマークの解析モデル・回転 Stokes 解との比較・誤差評価は [docs/sec6/iblbm2cdfSRT.md](../../docs/sec6/iblbm2cdfSRT.md)（Direct Forcing + SRT）、[docs/sec6/iblbm2cicMRT.md](../../docs/sec6/iblbm2cicMRT.md)（Implicit Correction + MRT）、[docs/sec6/iblbm2cicTRT.md](../../docs/sec6/iblbm2cicTRT.md)（Implicit Correction + TRT）を参照してください。接線速度の相対 L2 誤差は、DF-SRT 版が約 9.9%（$\tau=0.6$）、IC-MRT 版が約 7.7%（$\tau=0.6$）、IC-TRT 版が約 8.3%（$\tau_+=10$, $\Lambda=9/8$）です。陰的補正の 2 版は半力補正を含む Guo (2002) 力項により、Lagrangian 点が粗い（外 21／内 14 点）にもかかわらず DF-SRT より誤差が小さくなります。TRT 版は collision を偶奇 2 緩和に替え、磁数 $\Lambda$ で境界 slip を制御する点が特徴です。

### 粒子落下

| ファイル | 物理問題 |
| --- | --- |
| [iblbmsingle.c](iblbmsingle.c) | 単一円板（円柱断面）の重力沈降。終端速度と Stokes 解との比較 |
| [iblbmdkt.c](iblbmdkt.c) | 2 粒子の drafting–kissing–tumbling (DKT) ベンチマーク。後方粒子が前方粒子の wake に引き込まれて並走 → 接近 → 入れ替わる経典問題 |

> [iblbm2cdfSRT.c](iblbm2cdfSRT.c)、[iblbm2cicMRT.c](iblbm2cicMRT.c)、[iblbm2cicTRT.c](iblbm2cicTRT.c) には解説ドキュメントと可視化スクリプト ([scripts/plot_iblbm2cdfSRT_schematic.py](../../scripts/plot_iblbm2cdfSRT_schematic.py), [scripts/plot_iblbm2cdfSRT_results.py](../../scripts/plot_iblbm2cdfSRT_results.py), [scripts/plot_iblbm2cicMRT_schematic.py](../../scripts/plot_iblbm2cicMRT_schematic.py), [scripts/plot_iblbm2cicMRT_results.py](../../scripts/plot_iblbm2cicMRT_results.py), [scripts/plot_iblbm2cicTRT_schematic.py](../../scripts/plot_iblbm2cicTRT_schematic.py), [scripts/plot_iblbm2cicTRT_results.py](../../scripts/plot_iblbm2cicTRT_results.py)) を整備済みです。他のソースは冒頭ヘッダコメントに変数定義と離散化が記載されています。

## ビルドと実行

リポジトリのルートから次のコマンドで実行できます。出力先は既定で `outputs/sec6/<実行ファイル名>/` です。

```powershell
cmd /c scripts\run_one.cmd src\sec6\iblbm2cdfSRT.c
cmd /c scripts\run_one.cmd src\sec6\iblbmdkt.c
```

ビルドだけ行いたいときは [scripts/build_one.cmd](../../scripts/build_one.cmd) を使います。粒子落下系は格子が大きく、ヘッダコメントにもあるように「ワークステーション級の環境を前提」「Cygwin 等では格子を小さくする必要」との注記がついています。

## このディレクトリで扱う物理と数値

- **Immersed Boundary 法**：流体は構造化格子（Eulerian 格子）で解き、物体境界は Lagrangian 点列で別途表現します。各 Lagrangian 点で物体速度と流体補間速度の差から要する体積力 $\mathbf{F}_e$ を計算し、Dirac デルタの離散化（典型的には Roma–Peskin の 3-points hat 関数や 4-points cosine 関数）で Eulerian 格子へ分配します。
- **Direct Forcing (DF)**：補間 → 力計算 → 分配を陽的に行う最も単純な定式化。境界での速度残差が次ステップに残る。
- **Implicit Correction (IC)**：物体上の力分布を満たすように暗算的に反復補正を入れる、より精度の高い変種。
- **衝突演算子の選択**：SRT は単一緩和時間で最も基本。MRT は応力モーメントだけを物理粘性で緩和し、他のモーメントは安定性のために独立に調整。TRT は対称・非対称モーメントの 2 緩和時間（$\tau_+, \tau_-$）で「魔の数」関係 $\Lambda = (\tau_+ - 1/2)(\tau_- - 1/2)$ を一定に保つことで bounce-back 壁面位置の格子依存性を緩和。
- **粒子–流体–粒子相互作用**：[iblbmdkt.c](iblbmdkt.c) では粒子間の衝突応答（短距離反発 + 摩擦）も実装されており、複数粒子流の典型挙動である DKT を再現できます。

このセクションは、sec1–sec5 が固定境界条件・解析的境界形状を扱っていたのに対して、**任意形状・運動境界**を扱える IB-LBM へ拡張する位置づけです。
