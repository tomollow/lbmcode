# LBMcode

格子ボルツマン法に関する C のサンプルコードを、章ごとに整理したリポジトリです。

## 構成

- `src/sec1`: Taylor vortex の基本例
- `src/sec2`: FDM と LBM の追加例
- `src/sec3`: 自然対流と熱流動の例
- `src/sec4`: block と Cahn-Hilliard 関連の例
- `src/sec5`: Laplace と Zalesak の例
- `src/sec6`: immersed boundary method 関連の例

## ビルドと実行

現時点では、`lbmtv.c` 用のビルドスクリプトを同梱しています。

### Visual Studio を使う場合

リポジトリのルートで次を実行します。

```cmd
scripts\build_lbmtv.cmd
build\lbmtv.exe
```

このスクリプトは `vcvars64.bat` 経由で `cl.exe` を呼び出します。`lbmtv.c` ではスタック上に大きな配列を確保しているため、リンク時にスタックサイズも拡張しています。

### 出力ファイル

`build\lbmtv.exe` を実行すると、リポジトリのルートに次の結果ファイルを生成します。

- `error`
- `datautv`
- `datavtv`
- `datartv`
- `datautve`
- `datavtve`
- `datartve`

これらの生成物は Git の追跡対象から除外しています。

## 継続的インテグレーション

GitHub Actions では、`main` への push と pull request を契機に `src/sec1/lbmtv.c` をビルドします。