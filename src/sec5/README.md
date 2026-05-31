# sec5: 相場 LBM — 界面ダイナミクスと相場移流

本ディレクトリには、Cahn–Hilliard 型の相場（order parameter）方程式を D2Q9 LBM で解くサンプルが入っています。表面張力、界面厚さ、化学ポテンシャルを double-well ポテンシャルで表現し、流体運動（速度場）と相場（界面位置）を 2 つの分布関数で結合します。

## サンプル一覧

| ファイル | 物理問題 |
| --- | --- |
| [lbmlap.c](lbmlap.c) | Laplace の法則の検証：静止液滴内外の圧力差 $\Delta p = \sigma/R$ を半径をふってベンチマーク（解説: [docs/sec5/lbmlap.md](../../docs/sec5/lbmlap.md)） |
| [lbmzalesak.c](lbmzalesak.c) | Zalesak の円盤（切欠き付き円板の剛体回転）：相場移流の形状保持と数値拡散の評価（解説: [docs/sec5/lbmzalesak.md](../../docs/sec5/lbmzalesak.md)） |

> [lbmlap.c](lbmlap.c) は [docs/sec5/lbmlap.md](../../docs/sec5/lbmlap.md)、[lbmzalesak.c](lbmzalesak.c) は [docs/sec5/lbmzalesak.md](../../docs/sec5/lbmzalesak.md) に解説・図・ベンチマークを整備済みです。

## ビルドと実行

リポジトリのルートから次のコマンドで実行できます。出力先は既定で `outputs/sec5/<実行ファイル名>/` です。

```powershell
cmd /c scripts\run_one.cmd src\sec5\lbmlap.c
cmd /c scripts\run_one.cmd src\sec5\lbmzalesak.c
```

ビルドだけ行いたいときは [scripts/build_one.cmd](../../scripts/build_one.cmd) を使います。

## 可視化・ベンチマークスクリプト

次のスクリプトを用意しています（出力は `outputs/sec5/<実行ファイル名>/` と `docs/assets/sec5/` の両方）。

```powershell
# lbmlap.c (Laplace の法則)
d:/work/LBMcode/.venv/Scripts/python.exe scripts/plot_lbmlap_schematic.py          # 図 5.0 模式図
d:/work/LBMcode/.venv/Scripts/python.exe scripts/plot_lbmlap_results.py            # 図 5.1 結果 4 パネル
d:/work/LBMcode/.venv/Scripts/python.exe scripts/run_lbmlap_radius_sweep.py        # 図 5.2 半径スイープ

# lbmzalesak.c (Zalesak の円盤)
d:/work/LBMcode/.venv/Scripts/python.exe scripts/plot_lbmzalesak_schematic.py      # 図 5.0 模式図
d:/work/LBMcode/.venv/Scripts/python.exe scripts/plot_lbmzalesak_results.py        # 図 5.1 1 回転移流 4 パネル
d:/work/LBMcode/.venv/Scripts/python.exe scripts/run_lbmzalesak_peclet_sweep.py    # 図 5.2 Peclet スイープ
```

## このディレクトリで扱う物理と数値

- **相場法（Cahn–Hilliard 型）**：界面を厚さ $W$ の遷移層として表現する diffuse-interface 法。化学ポテンシャル
  $$\mu = 4\,\beta\,\phi\,(\phi^2 - \phi_0^2) - \kappa\,\nabla^2\phi$$
  と移動度 $M$ により相場の輸送 $\partial_t\phi + \mathbf{u}\cdot\nabla\phi = M\nabla^2\mu$ を解きます。表面張力 $\sigma$ は、本コードでは内部で $\beta = (3/4)\sigma\phi_0^4/W$, $\kappa = (3/8)\sigma W/\phi_0^2$ と入力 $\sigma$ から逆算しており、$\sigma$ を直接指定できる構成になっています。
- **双分布関数構造**：lbmnc.c などと同じく、速度場（D2Q9, 分布関数 $f$）と相場（D2Q9, 分布関数 $g$）を別個に解き、毎ステップ巨視量を再構成します。表面張力は化学ポテンシャル勾配からの体積力 $\mathbf{F} = \mu\nabla\phi$ として速度分布に加えます。
- **Laplace の法則**：半径 $R$ の液滴に対し圧力差 $\Delta p = \sigma/R$ が成り立つかを、$R$ や格子解像度を変えて検証する古典ベンチマーク（[lbmlap.c](lbmlap.c)）。
- **Zalesak 円盤**：速度場を解かず解析的な剛体回転場を与え、切欠き円板を**ちょうど 1 回転**（2500 ステップ）させて初期形状との一致を見る古典的移流テスト（[lbmzalesak.c](lbmzalesak.c)）。相場移流スキーム単体の数値拡散・分散・保存性を、幾何誤差 $E_1$ や Peclet 数依存性として評価できます。

このセクションは sec3 の D2Q5 熱輸送ベンチマーク（lbmtherm.c）と同様、解析解や閉形式ベンチマークと直接比較できる検証問題群です。
