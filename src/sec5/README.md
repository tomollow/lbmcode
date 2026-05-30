# sec5: 相場 LBM — 界面ダイナミクスと相場移流

本ディレクトリには、Cahn–Hilliard 型の相場（order parameter）方程式を D2Q9 LBM で解くサンプルが入っています。表面張力、界面厚さ、化学ポテンシャルを double-well ポテンシャルで表現し、流体運動（速度場）と相場（界面位置）を 2 つの分布関数で結合します。

## サンプル一覧

| ファイル | 物理問題 |
| --- | --- |
| [lbmlap.c](lbmlap.c) | Laplace の法則の検証：静止液滴内外の圧力差 $\Delta p = \sigma/R$ を界面厚さと格子解像度をふってベンチマーク |
| [lbmzalesak.c](lbmzalesak.c) | Zalesak の円盤（切欠き付き円板の回転）：相場移流の形状保持と数値拡散の評価 |

> 本セクションには `docs/sec5/` の解説ドキュメントおよび専用の可視化スクリプトはまだ整備されていません。各ソースの冒頭ヘッダコメントに変数の意味と離散化が記載されているので、まずはそちらを参照してください。

## ビルドと実行

リポジトリのルートから次のコマンドで実行できます。出力先は既定で `outputs/sec5/<実行ファイル名>/` です。

```powershell
cmd /c scripts\run_one.cmd src\sec5\lbmlap.c
cmd /c scripts\run_one.cmd src\sec5\lbmzalesak.c
```

ビルドだけ行いたいときは [scripts/build_one.cmd](../../scripts/build_one.cmd) を使います。

## このディレクトリで扱う物理と数値

- **相場法（Cahn–Hilliard 型）**：界面を厚さ $W$ の遷移層として表現する diffuse-interface 法。化学ポテンシャル
  $$\mu = 4\,\beta\,\phi\,(\phi^2 - \phi_0^2) - \kappa\,\nabla^2\phi$$
  と移動度 $M$ により相場の輸送 $\partial_t\phi + \mathbf{u}\cdot\nabla\phi = M\nabla^2\mu$ を解きます。表面張力 $\sigma$ は、本コードでは内部で $\beta = (3/4)\sigma\phi_0^4/W$, $\kappa = (3/8)\sigma W/\phi_0^2$ と入力 $\sigma$ から逆算しており、$\sigma$ を直接指定できる構成になっています。
- **双分布関数構造**：lbmnc.c などと同じく、速度場（D2Q9, 分布関数 $f$）と相場（D2Q9, 分布関数 $g$）を別個に解き、毎ステップ巨視量を再構成します。表面張力は化学ポテンシャル勾配からの体積力 $\mathbf{F} = \mu\nabla\phi$ として速度分布に加えます。
- **Laplace の法則**：半径 $R$ の液滴に対し圧力差 $\Delta p = \sigma/R$ が成り立つかを、$R$ や格子解像度を変えて検証する古典ベンチマーク（[lbmlap.c](lbmlap.c)）。
- **Zalesak 円盤**：剛体回転場で形状を保ったまま回り続けるべき切欠き円板の数値解像を見るテスト（[lbmzalesak.c](lbmzalesak.c)）。格子粗化や移動度の選び方による界面のぼやけと体積保存性を評価できます。

このセクションは sec3 の D2Q5 熱輸送ベンチマーク（lbmtherm.c）と同様、解析解や閉形式ベンチマークと直接比較できる検証問題群です。
