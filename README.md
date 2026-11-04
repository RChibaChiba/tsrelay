# tsrelay

基質（生成物）と触媒の xyz、settings.ini だけを入力に、xtb dock から UMA による TS の振動解析までを自走させる。
人の操作は「起動」と「dock 結果の選択」の 2 回だけ。

```
WSL（ローカル）                               Colab（無料 T4）
入力検証 → コア検出 → xtb dock → 選択 GUI  →  UMA で p 錯体 → 逆スキャン → sm 錯体 → DMF → Sella → 振動解析
```

## 前提

- Windows の WSL2（Ubuntu）
- xtb 6.6.1（6.7.1 の公式バイナリは `xtb dock` が segfault する）
- [google-colab-cli](https://pypi.org/project/google-colab-cli/)（`colab --auth oauth2 sessions` が通ること）
- Hugging Face のトークン（facebook/UMA の利用許可を得たアカウントで `hf auth login`、または環境変数 `HF_TOKEN`）
- UMA の重みは手元に要らない。ランタイムが HF から直接取る（手元から 2 GB 級を送ると遅いため）

## 使い方

```bash
uv venv .venv && uv pip install --python .venv/bin/python ase numpy scipy
cd examples/pd_re
../../.venv/bin/python ../../run.py settings.ini
```

1. 起動すると、検出したコアとスキャンするペアを表示して y/n を聞く（`--yes` で省略）。
2. xtb dock の候補がブラウザに並ぶ（コアは灰色、directed 原子は赤）。選んで「確定」を押す。
3. T4 のランタイムを作り、環境の準備・重みの取得・doctor を済ませてから、選んだ構造ごとに TS を求める。
   済んだものから `out/<ID>/` に回収し、全部済むと `out/summary.csv` を書いてランタイムを止める。

| オプション | 意味 |
| --- | --- |
| `--yes` | 起動時の y/n 確認を省く |
| `--keep` | 終わってもランタイムを止めない（調べるとき用。止めるときは `colab stop -s <session>`） |
| `--dock XYZ` | 手元の dock の結果を使い、xtb dock を飛ばす（xtb がなくても動く） |

`--dock` には `xtb dock` の `final_structures.xyz` を渡す。

```bash
../../.venv/bin/python ../../run.py settings.ini --dock final_structures.xyz
```

- 同じ形式なら別の手段で作った構造でもよい。原子は「基質が先、触媒が後」の順で、コメント行の先頭はエネルギー（Hartree）。
- 候補の検査は通常と同じ。コアが settings.ini の触媒から 0.05 Å 以上ずれた構造は捨てるので、コアを固定して作った構造を渡す。

- 途中で止まっても、同じコマンドで続きから再開する。段階をやり直すときは `work/NN_*` を消す。
- T4 が割り当てられないときは 30 分おきに 3 回試し、だめなら理由を表示して終了する。
- settings.ini を変えて同じ workdir で回そうとすると、差分を表示して止まる。

## settings.ini

番号はすべて 1 始まり。錯体の番号は「基質が先、触媒が後」。例は [examples/pd_re/settings.ini](examples/pd_re/settings.ini)。

| 節 | 主な項目 |
| --- | --- |
| `[input]` | `substrate`、`catalyst`、`workdir` |
| `[state]` | 全体・基質・触媒の電荷と多重度 |
| `[dock]` | `directed_atoms`（金属に向ける基質原子）、`xtb_args` |
| `[core]` | `cutoff_scale`、`add_atoms` / `drop_atoms`（触媒の番号）。金属が複数あればすべてコアにする |
| `[scan]` | `bonds = i-j:目標Å, ...`（伸ばすペアと縮めるペア）、`nstep` |
| `[uma]` | `model`（facebook/UMA のモデル名。既定は `uma-s-1p2p1`） |
| `[ts]` | `fmax_opt`、`fmax_ts` |
| `[run]` | `cores`（xtb のスレッド数）、`session`（Colab のランタイム名） |

## 出力

```
workdir/
  settings.ini            起動時の設定の写し
  log/tsrelay.log
  work/00_check 01_dock 02_pick 03_remote   段階ごとの作業場所（done で完了）
  out/summary.csv         id, status, E_sm_eV, E_ts_eV, E_p_eV, dE_kcal, imag_cm（ok を先に、ΔE‡ の小さい順）
  out/<ID>/               p.xyz scan.xyz sm.xyz path.xyz ts_guess.xyz ts.xyz ts_mode.xyz freq.txt sella.log result.json
```

status が ok になるのは、Sella が収束し、50 cm⁻¹ 以上の虚数振動がちょうど 1 つのとき。
虚数モードがスキャンした原子をほとんど動かさない場合などは、warnings に注記する。

## 例（examples/pd_re）

Pd(PH₃)₂ からのエタンの還元的脱離を、生成物側（エタン錯体）から逆にたどる 17 原子の系。
T4 で ΔE‡ ≈ 26 kcal/mol、虚数振動 −478 cm⁻¹ が 1 つ。1 構造 30 秒〜2 分程度、
ランタイムごとに環境の準備に約 4 分、重みの取得と doctor に約 3 分かかる
（`--keep` で残したランタイムを再利用すれば、準備と取得は省かれる）。

## ファイル

| ファイル | 実行場所 | 役割 |
| --- | --- | --- |
| run.py | WSL | 入口。段階を順に呼び、done で飛ばす。結果の回収と summary |
| local.py | WSL | ini の読込と検証、コア検出、xtb dock、候補の後処理 |
| pick.py | WSL | 選択 GUI（127.0.0.1 の小さな HTTP サーバー） |
| colab.py | WSL | colab CLI の薄いラッパー（目印行の読み取り、再試行、障害の判定） |
| remote.py | Colab | 構造ごとの UMA 計算、doctor、重みの取得（`python remote.py bundle_dir` で単体でも動く） |

仕様の詳細は「tsrelay 仕様書」（Claude Docs）を参照。

## ライセンス

このリポジトリのコードは [MIT](LICENSE)。UMA の重みは含まず、facebook/UMA の利用規約に従って各自で取得する。
xtb、fairchem-core、pydmf、Sella などの依存先はそれぞれのライセンスに従う。
