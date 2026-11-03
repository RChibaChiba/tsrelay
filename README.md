# tsrelay

基質（生成物）と触媒の xyz、settings.ini だけを入力に、xtb dock から UMA による TS の振動解析までを自走させる。
人の操作は「起動」と「dock 結果の選択」の 2 回だけ。

```
WSL（ローカル 6 コア）                         Colab（無料 T4）
入力検証 → コア検出 → xtb dock → 選択 GUI  →  UMA で p 錯体 → 逆スキャン → sm 錯体 → DMF → Sella → 振動解析
```

## 前提

- Windows の WSL2（Ubuntu）
- xtb 6.6.1（6.7.1 の公式バイナリは `xtb dock` が segfault する）
- [google-colab-cli](https://pypi.org/project/google-colab-cli/)（`colab --auth oauth2 sessions` が通ること）
- UMA の重み `uma-s-1p2p1.pt`（Hugging Face の facebook/UMA から手元に一度ダウンロードしておく。ランタイムでは HF にログインしない）

## 使い方

```bash
uv venv .venv && uv pip install --python .venv/bin/python ase numpy scipy
cd examples/pd_re
../../.venv/bin/python ../../run.py settings.ini
```

1. 起動すると、検出したコアとスキャンするペアを表示して y/n を聞く（`--yes` で省略）。
2. xtb dock の候補がブラウザに並ぶ。選んで「確定」を押す。
3. 選んだ構造ごとに Colab で TS を求め、済んだものから `out/<ID>/` に回収する。全部済むと `out/summary.csv` を書く。

途中で止まっても、同じコマンドで続きから再開する。段階をやり直すときは `work/NN_*` を消す。

## ファイル

| ファイル | 実行場所 | 役割 |
| --- | --- | --- |
| run.py | WSL | 入口。段階を順に呼び、done で飛ばす |
| local.py | WSL | ini の読込と検証、コア検出、xtb dock、候補の後処理 |
| pick.py | WSL | 選択 GUI（127.0.0.1 の小さな HTTP サーバー） |
| colab.py | WSL | colab CLI の薄いラッパー（大きなファイルの分割送付、障害の判定） |
| remote.py | Colab | 構造ごとの UMA 計算（`python remote.py bundle_dir` で単体でも動く） |

仕様の詳細は「tsrelay 仕様書」（Claude Docs）を参照。
