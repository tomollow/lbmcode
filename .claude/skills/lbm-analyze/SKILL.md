---
description: Document, run, analyse, and benchmark an LBM C sample in src/secN/ — produces docs/secN/<name>.md, schematic and result figures, benchmark CSVs, parameter-sweep scripts, and a PR. Use when the user wants to "analyse / document / benchmark" any C source under src/secN/, or says "解説して", "解析して", "ベンチマーク比較して" against any sample in this repo.
---

# LBM サンプル解説・解析パイプライン

このリポジトリ (lbmcode) の `src/secN/*.c` を、ドキュメント・図・ベンチマーク比較・PR まで一式作成する手順です。

## 入力の確認

ユーザーから以下を確認 (不明なら聞く)：

- **対象ファイル**: `src/sec{N}/{name}.c`
- **比較するベンチマーク**: 文献名 (例: de Vahl Davis 1983, Ghia 1982) または「コード内の解析解」
- **すでにある成果物**: `docs/sec{N}/{name}.md` や `scripts/plot_{name}*.py` が既にあれば、上書きでなく追補する

## 成果物

1. `docs/sec{N}/{name}.md` — 日本語の解説ドキュメント
2. `docs/assets/sec{N}/{name}_schematic.png` — 解析モデル模式図
3. `docs/assets/sec{N}/{name}_results.png` — 主要結果図 (4–5 panel)
4. `docs/sec{N}/generated/{name}_*.csv` — 数値ベンチマーク比較表
5. `scripts/plot_{name}_schematic.py` — 模式図再生成
6. `scripts/plot_{name}_results.py` — 結果図再生成
7. (任意) `scripts/run_{name}_*_sweep.py` — パラメータスイープ
8. `src/sec{N}/README.md` — セクション全体オリエンテーション (存在しなければ新規)
9. PR (feature branch `claude/sec{N}-<topic>`)

## ワークフロー (順番厳守)

### Phase 1 — 解説ドキュメント

`docs/sec2/lbmcavi.md` を雛形に同じ構成で書く:

```
概要 → 扱う物理量 (記号表) → 解析モデル (図 N.0 を参照)
→ 対象問題 (支配方程式) → 格子モデル (D2Q9/D2Q5 の重み・離散速度)
→ 平衡分布関数 → 時間発展 (collision → forcing → streaming → 巨視量再構成)
→ 境界条件 → 数値設定 → 解析結果 → このコードの見どころ
```

- 数式は LaTeX
- ソース該当行に `[name.c:line](../../src/sec{N}/{name}.c#L<line>)` リンク
- 主要な定数 (重み, 緩和率, 平衡分布係数, 浮力係数) は **コード値と数式が一致することを逐行確認**

### Phase 2 — 解析モデル図

`scripts/plot_{name}_schematic.py` を作成し、`docs/assets/sec{N}/{name}_schematic.png` を生成。

含めるもの: 計算領域, 境界条件 (色/ハッチで識別), 外力方向, 座標系, 代表長さ寸法線, 期待される典型流動パターン (矢印で), 主要無次元数 (Ra, Pr, Re 等)。

**スクリプト共通テンプレ**:

```python
ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec{N}" / "{name}"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec{N}"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                    "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
```

出力は **両方** に保存: `outputs/sec{N}/{name}/{filename}.png` と `docs/assets/sec{N}/{filename}.png`。

### Phase 3 — 実行 + 結果プロット

1. ビルド + 実行:
   ```powershell
   cmd /c scripts\run_one.cmd src\sec{N}\{name}.c
   ```
2. `outputs/sec{N}/{name}/` の data ファイルを numpy で読み込み
3. 4–5 パネル図に整理:
   - フィールド (T, ω, ψ 等の contour)
   - 中心線プロファイル (ベンチマーク値を ★ で重ね描き)
   - 局所診断量 (Nu(y), C_d(t), スペクトル 等)
   - 誤差ヒストグラム (解析解と比較する場合)
4. **ベンチマーク比較表** を CSV と Markdown table の両方で出力
5. 一致度を %誤差で正直に報告 (誇張禁止)

### Phase 4 — Next Action ループ (順番厳守)

ベンチマークとの差が大きい量があれば、以下の順で改善:

1. **出力解像度フル化** — C ソースの出力ループの stride を 1 にし、物理内部 `i,j = 1..n-1` を全点書き出す (`fopen("data...")` 部分のみ修正、**計算ロジックは触らない**)
2. **数値評価の高精度化** — 例: 壁面勾配を 2 点 → 3 点片側 2 次差分。半セルずれの half-way bounce-back では:
   ```
   dT/dx|_wall = (-8 T_w + 9 T_1 - T_2) / (3 h)
   ```
   ($T_w$ は壁温, $T_1, T_2$ は壁から数えて 1, 2 番目の内部セル中心, $h = 1/n$ は無次元セル幅)
3. **パラメータスイープ** — 制御パラメータを 4 値程度で sweep し、ベンチマーク表と log-log 比較。発散ケースは `"diverged": 1.0` として除外し、CSV/プロットでスキップ。一時 `.c` と `.exe` は `try/finally` で必ず削除。`subprocess.run` には `timeout` を入れる
4. **兄弟ファイルへの展開** — 同セクションに関連ファイルがあれば doc + section README に追加

### Phase 5 — 自己レビュー (必須)

`Agent(subagent_type="general-purpose", ...)` で別エージェントに critical review を依頼:

```
PR の独立クリティカルレビュー。以下を疑え:
- 数式の係数を再導出して一致するか
- 物理直感 (回転方向, ピーク位置) と数値結果が整合するか
- argmax / argmin の符号付き/絶対値の規約はベンチマーク文献と一致するか
  (DVD/Ghia 規約は signed argmax の正のピーク。np.abs で tie-break しない)
- regex 置換が想定外マッチしないか (行頭アンカ ^\s* または \b で固定)
- 例外時に一時ファイル/exe が orphan しないか (try/finally で守る)
- ドキュメント記述と現在のコードが食い違っていないか
  (特に修正前の古い文章が残っていないか)
- 回転方向: 左 hot / 右 cold + 重力 -y → 時計回り (CW)。
  これは間違えやすいので最初に符号テストする
```

報告された致命的問題は必ず修正してから commit。軽微なものは fixup commit で対応可。

### Phase 6 — Commit, PR, レビュー対応, マージ

1. ブランチ作成:
   ```bash
   git checkout -b claude/sec{N}-<topic>
   ```
2. ステージ。**`.claude/` は手動で除外** (skill 以外は local state):
   ```bash
   git add src/sec{N}/{name}.c src/sec{N}/README.md \
           docs/sec{N}/ docs/assets/sec{N}/ \
           scripts/plot_{name}_*.py scripts/run_{name}_*.py
   ```
3. コミット (Co-Authored-By 必須、HEREDOC で改行を含む):
   ```bash
   git commit -m "$(cat <<'EOF'
   Add sec{N} docs: {name}.c <one-line description>

   <Multi-line body describing what changed, why, what numbers reproduce
   the benchmark, and any known limitations>

   Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
   EOF
   )"
   ```
4. push & PR:
   ```bash
   git push -u origin claude/sec{N}-<topic>
   gh pr create --title "..." --body "..."
   ```
   PR body は日本語、`## Summary` と `## Test plan` のチェックリスト形式。
5. レビュー指摘は fixup commit で。マージ可否は `gh pr view <N> --json mergeable,mergeStateStatus,statusCheckRollup` で確認 (CLEAN かつ build SUCCESS)
6. マージ (過去 PR の慣習に合わせて squash でも rebase でもなく merge commit):
   ```bash
   gh pr merge <N> --merge
   ```
7. `git checkout main && git pull --ff-only`

## 規約

| 項目 | 規約 |
| --- | --- |
| 言語 | docs は日本語、コード/コメント/コミット/PR タイトルは英語 |
| 図番号 | 図 {N}.0 = 解析モデル, {N}.1 = 主要結果, {N}.2 以降 = 追加診断 |
| スクリプト出力先 | `outputs/sec{N}/{name}/` と `docs/assets/sec{N}/` の両方 |
| CSV 配置 | `docs/sec{N}/generated/` |
| section README | sec3 の README.md を雛形に: ファイル一覧表 + ビルド/実行 + スクリプト一覧 + 「このディレクトリで扱う物理」 |
| ベンチマーク規約 | $u_{\max}$ は鉛直中心線の **正のピーク** (`np.argmax(u_vert)`), $v_{\max}$ は水平中心線の **正のピーク**。`np.abs` で tie-break させない |
| `subprocess.run` | 必ず `shell=False, timeout=<sec>`、`shell=True` 禁止 |
| regex 置換 | `(?m)^\s*…` で行頭アンカ、`count=1` を併用 |

## 起動例

```
/lbm-analyze
対象: src/sec5/lbmlap.c
ベンチマーク: Laplace の法則 Δp = σ/R
特記事項: 界面厚さ W と σ/(σ_theory) の依存性も評価したい
```

## チェックリスト (PR 提出前)

- [ ] 全数式の係数をコードと逐行確認した
- [ ] ベンチマーク数値の出典を URL or 文献名で明記した
- [ ] 図のキャプションが図番号付きで、本文から参照されている
- [ ] スクリプトを clean state で再実行し、生成物が再現する
- [ ] regex / 一時ファイル削除が `try/finally` で例外安全
- [ ] subprocess に `timeout` 設定
- [ ] 回転方向・ピーク位置などの物理直感を符号テストで確認した
- [ ] サブエージェントの critical review を一度通した
