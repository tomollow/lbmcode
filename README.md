# LBMcode

格子ボルツマン法に関する C のサンプルコードを、章ごとに整理したリポジトリです。

## 構成

- `src/sec1`: Taylor vortex の基本例
- `src/sec2`: Poiseuille flow、cavity flow、移流方程式などの基本例
- `src/sec3`: 熱 LBM と自然対流の例
- `src/sec4`: multi-block と cavity flow の例
- `src/sec5`: 二相流に関する Laplace 則と Zalesak disk の例
- `src/sec6`: immersed boundary method を使った円筒 Couette flow などの例
- `docs/sec1`, `docs/sec2`: 各サンプルコードの説明ドキュメント
- `docs/assets`: 図生成スクリプトが出力する公開用画像
- `docs/plot_generation.md`: 図生成スクリプト一覧と再生成方法

### ファイル一覧

- `src/sec1/lbmtv.c`: Taylor vortex flow
- `docs/sec1/lbmtv.md`: `lbmtv.c` の説明と数式
- `src/sec2/fdlbm.c`: 有限差分格子ボルツマン法による Poiseuille flow
- `docs/sec2/fdlbm.md`: `fdlbm.c` の説明と数式
- `src/sec2/fdmadv.c`: 1 次元移流方程式の有限差分法
- `docs/sec2/fdmadv.md`: `fdmadv.c` の説明と数式
- `src/sec2/lbmbound.c`: 境界条件スキームを比較する Poiseuille flow
- `docs/sec2/lbmbound.md`: `lbmbound.c` の説明と数式
- `src/sec2/lbmcavi.c`: compressibility error を含む cavity flow
- `docs/sec2/lbmcavi.md`: `lbmcavi.c` の説明と数式
- `src/sec2/lbmpoi.c`: 格子ボルツマン法による Poiseuille flow
- `docs/sec2/lbmpoi.md`: `lbmpoi.c` の説明と数式
- `src/sec3/lbmnc.c`: double-population thermal LBM による自然対流
- `src/sec3/lbmtherm.c`: Thermal LBM の基本例
- `src/sec4/lbmblock.c`: multi-block LBM による Couette flow
- `src/sec4/lbmcm.c`: cavity flow の例
- `src/sec5/lbmlap.c`: Laplace's law の検証
- `src/sec5/lbmzalesak.c`: Zalesak's disk の移流
- `src/sec6/iblbm2cdfSRT.c`: direct forcing 法による円筒 Couette flow
- `src/sec6/iblbm2cicMRT.c`: implicit correction 法と MRT による円筒 Couette flow
- `src/sec6/iblbm2cicTRT.c`: implicit correction 法と TRT による円筒 Couette flow
- `src/sec6/iblbmdkt.c`: immersed boundary-lattice Boltzmann method の例
- `src/sec6/iblbmsingle.c`: immersed boundary-lattice Boltzmann method の単体例

## ビルドと実行

Visual Studio の C/C++ ツールチェーンを使う前提で、単体ビルドと一括ビルドのスクリプトを同梱しています。

### Visual Studio を使う場合

リポジトリのルートで次を実行します。

```cmd
scripts\build_lbmtv.cmd
scripts\run_one.cmd src\sec1\lbmtv.c
```

任意の 1 ファイルだけビルドする場合は次を使えます。

```cmd
scripts\build_one.cmd src\sec2\fdmadv.c
build\bin\fdmadv.exe
```

実行生成物をルートに散らさずに実行する場合は、次を使えます。

```cmd
scripts\run_one.cmd src\sec2\fdmadv.c
```

この場合、実行ファイルは `build\bin` に、実行時に生成されるファイルは `outputs\sec2\fdmadv` のように sec ごとのフォルダに出力されます。

すべての C ファイルをまとめてビルドする場合は次を実行します。

```cmd
scripts\build_all.cmd
```

生成物をまとめて消して作業ディレクトリをきれいに戻す場合は、次を実行します。

```cmd
scripts\clean.cmd
```

これらのスクリプトは必要に応じて `vcvars64.bat` を呼び出し、`cl.exe` を使ってビルドします。いくつかのコードはスタック上に大きな配列を確保しているため、リンク時にスタックサイズも拡張しています。

### 出力ファイル

`scripts\run_one.cmd src\sec1\lbmtv.c` を実行すると、`outputs\sec1\lbmtv` に次の結果ファイルを生成します。

- `error`
- `datautv`
- `datavtv`
- `datartv`
- `datautve`
- `datavtve`
- `datartve`

これらの生成物は Git の追跡対象から除外しています。

## 実行例

各 sec の代表例を試す最短コマンドです。

### sec1

```cmd
scripts\run_one.cmd src\sec1\lbmtv.c
```

### sec2

```cmd
scripts\run_one.cmd src\sec2\fdlbm.c
```

`fdmadv.c` は対話入力でもコマンドライン引数でも実行できます。たとえば Upwind scheme を既定の `dt = 0.5` で実行するなら次のどちらでも動きます。

```cmd
echo 1 | build\bin\fdmadv.exe
```

```cmd
build\bin\fdmadv.exe 1 0.5
```

この実行では最終プロファイルを `fdmadv` に、各表示時刻の履歴を `fdmadv_history.dat` に出力します。

Courant 数 `0.5` と `1.0` の比較図を再生成する場合は次を実行します。

```powershell
d:/work/LBMcode/.venv/Scripts/python.exe scripts/plot_fdmadv_courant_compare.py
```

生成される図は `outputs/sec2/fdmadv_courant_compare.png` に保存され、公開用コピーが `docs/assets/sec2/fdmadv_courant_compare.png` に作成されます。

図生成スクリプト全体の一覧は `docs/plot_generation.md` にまとめています。

### sec3

```cmd
scripts\run_one.cmd src\sec3\lbmtherm.c
```

### sec4

```cmd
scripts\run_one.cmd src\sec4\lbmblock.c
```

### sec5

```cmd
scripts\run_one.cmd src\sec5\lbmlap.c
```

### sec6

```cmd
scripts\run_one.cmd src\sec6\iblbm2cdfSRT.c
```

## 継続的インテグレーション

GitHub Actions では、`main` への push と pull request を契機に次を実行します。

- すべての C ファイルのビルド
- `lbmtv.exe` の smoke test 実行

## 更新履歴

`main` ブランチに取り込まれた主な変更（新しい順）：

- **2026-05-23** [#22](https://github.com/tomollow/lbmcode/pull/22): sec4 の [`src/sec4/lbmblock.c`](src/sec4/lbmblock.c) に 2D 場ダンプ（t=3000 で粗・細格子それぞれの $u, v$）を追加し、多重格子可視化スクリプト [`scripts/plot_lbmblock_couette.py`](scripts/plot_lbmblock_couette.py) と [`scripts/plot_lbmblock_distribution.py`](scripts/plot_lbmblock_distribution.py) を同梱。[`docs/sec4/lbmblock.md`](docs/sec4/lbmblock.md) にメッシュ配置図・分布結果（$t=3000$）・過渡応答セクションを追加し、Filippova-Hänel スケーリング則とインターフェース処理が解析解（最大相対誤差 $\approx 7\times 10^{-3} u_t$）まで収束することを確認
- **2026-05-20** [#21](https://github.com/tomollow/lbmcode/pull/21): sec4 の [`src/sec4/lbmcm.c`](src/sec4/lbmcm.c) 説明ドキュメント [`docs/sec4/lbmcm.md`](docs/sec4/lbmcm.md) を追加。SRT / MRT / 中心モーメント (CM) 衝突の数式（モーメント変換 $M$、シフト行列 $N(\mathbf{u})$、緩和行列 $S$）と Re=100/1000/5000 の比較結果を整理。6 ケースを一括実行するヘルパー [`scripts/run_lbmcm_compare.ps1`](scripts/run_lbmcm_compare.ps1) と分布・比較プロット [`scripts/plot_lbmcm_distribution.py`](scripts/plot_lbmcm_distribution.py) / [`scripts/plot_lbmcm_compare.py`](scripts/plot_lbmcm_compare.py) も同梱
- **2026-05-19** [#19](https://github.com/tomollow/lbmcode/pull/19): sec4 に Spalart-Allmaras DES97 のハイブリッド RANS/LES 実装を追加（cavity / karman / backward_step）。共通定数とクロージャ関数を [`src/sec4/sa_closure.h`](src/sec4/sa_closure.h) に集約。高 Re サニティ版 `karman_des_hires.c` と RANS/LES 領域マップ [`scripts/plot_des_region_map.py`](scripts/plot_des_region_map.py) も同梱。クロスケース比較 doc: [`docs/sec4/les_summary.md`](docs/sec4/les_summary.md), [`docs/sec4/keps_summary.md`](docs/sec4/keps_summary.md)
- **2026-05-12** [#18](https://github.com/tomollow/lbmcode/pull/18): sec4 の 6 ケース（kelbm / taylor_green / kelvin_helmholtz / cavity / backward_step / karman）に Smagorinsky LES 版を追加し、pure / k-ε / LES の 3 variant 比較に再編。[`docs/sec4/les_summary.md`](docs/sec4/les_summary.md) でクロスケース比較を新設
- **2026-05-11** [#17](https://github.com/tomollow/lbmcode/pull/17): 同 6 ケースを横断する k-ε モデルの体系的比較 [`docs/sec4/keps_summary.md`](docs/sec4/keps_summary.md) を追加
- **2026-05-11** [#16](https://github.com/tomollow/lbmcode/pull/16): Karman 渦列ケース（pure LBM + k-ε）を追加（[`src/sec4/karman.c`](src/sec4/karman.c), [`src/sec4/karman_keps.c`](src/sec4/karman_keps.c)）
- **2026-05-11** [#15](https://github.com/tomollow/lbmcode/pull/15): cavity のコーナーセル壁関数を修正（4 壁の `apply_wall_function` で `max` ロジック導入）
- **2026-05-11** [#14](https://github.com/tomollow/lbmcode/pull/14): 後方ステップ流れケース（pure LBM + k-ε）を追加（[`src/sec4/backward_step.c`](src/sec4/backward_step.c), [`src/sec4/backward_step_keps.c`](src/sec4/backward_step_keps.c)）
- **2026-05-10** [#13](https://github.com/tomollow/lbmcode/pull/13): 蓋駆動 cavity ケース（pure LBM + k-ε）を追加（[`src/sec4/cavity.c`](src/sec4/cavity.c), [`src/sec4/cavity_keps.c`](src/sec4/cavity_keps.c)）
- **2026-05-10** [#12](https://github.com/tomollow/lbmcode/pull/12): Kelvin-Helmholtz 不安定ケース（pure LBM + k-ε）を追加
- **2026-05-10** [#11](https://github.com/tomollow/lbmcode/pull/11): Taylor-Green 渦ケース（pure LBM + k-ε）を追加
- **2026-05-10** [#10](https://github.com/tomollow/lbmcode/pull/10): チャンネル流の k-ε 乱流ケース [`src/sec4/kelbm.c`](src/sec4/kelbm.c) を追加