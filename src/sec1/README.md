# sec1: LBM の基礎 — Taylor 渦と解析解との比較

本ディレクトリには、格子ボルツマン法（LBM）の最も基本的な検証として、解析解が閉形式で得られる Taylor–Green 渦（Taylor vortex）を BGK 衝突演算子で解くプログラムが入っています。LBM の収束次数や格子解像度依存性を、解析解と直接比較して評価するための入口に位置づけられます。

## サンプル一覧

| ファイル | 目的 | ドキュメント |
| --- | --- | --- |
| [lbmtv.c](lbmtv.c) | 2 次元 Taylor 渦の D2Q9-BGK 数値解と粘性減衰する解析解との誤差評価 | [docs/sec1/lbmtv.md](../../docs/sec1/lbmtv.md) |

## ビルドと実行

リポジトリのルートから次のコマンドで実行できます。出力先は既定で `outputs/sec1/<実行ファイル名>/` です。

```powershell
cmd /c scripts\run_one.cmd src\sec1\lbmtv.c
```

ビルドだけ行いたいときは [scripts/build_one.cmd](../../scripts/build_one.cmd)、まとめてビルドする場合は [scripts/build_all.cmd](../../scripts/build_all.cmd) を使います。

## 解析・可視化スクリプト

| スクリプト | 出力 |
| --- | --- |
| [scripts/plot_lbmtv_nx_study.py](../../scripts/plot_lbmtv_nx_study.py) | 格子解像度 $n_x$ を変えたときの誤差収束次数の確認図 |

## このディレクトリで扱う物理

- **Taylor 渦の解析解**：粘性流体中で初期に与えた 2 重周期的な渦は、運動エネルギーを $e^{-2\nu k^2 t}$ で指数減衰させながら形状を保ちます。LBM の誤差源（離散化・圧縮性）と分けて評価しやすい代表問題です。
- **D2Q9-BGK**：9 速度・単一緩和時間モデルの最小構成。緩和時間 $\tau$ と動粘性係数 $\nu = (\tau - 0.5)/3$ の関係、平衡分布関数 $f_k^{\mathrm{eq}} = w_k \rho (1 + 3\mathbf{c}_k\cdot\mathbf{u} + \frac{9}{2}(\mathbf{c}_k\cdot\mathbf{u})^2 - \frac{3}{2}\lvert\mathbf{u}\rvert^2)$ の実装を、後続セクションのベースとして確認できます。
- **格子収束次数**：$n_x$ をスイープして $L_2$ 誤差を測ることで、空間 2 次精度であることを確認できます。

このセクションが後続セクション全体の起点で、ここで導入される BGK・D2Q9 の実装と緩和時間の取り扱いが、sec2 以降の境界条件、応用流れ、熱輸送、乱流モデルへと拡張されていきます。
