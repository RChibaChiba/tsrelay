"""tsrelay の入口。段階を順に呼び、done のある段階は飛ばす。

    python run.py settings.ini [--yes] [--keep]

段階ごとのフォルダ（workdir/work/NN_name）に done を置くのが唯一の状態管理。消せばその段階からやり直す。
"""
import argparse
import csv
import difflib
import json
import os
import shutil
import sys
import tarfile
import time
from datetime import datetime
from pathlib import Path

import colab
import local
import pick
from remote import PINS

HERE = Path(__file__).resolve().parent
POLL_S = 60
SETUP = ("if [ ! -f setup.ok ]; then apt-get -qq update && apt-get -qq install -y coinor-libipopt-dev pkg-config"
         " && pip install -q " + " ".join(f"{k}=={v}" for k, v in PINS.items()) + " && touch setup.ok; fi")
# トークンのファイルは取得の成否にかかわらず消す
FETCH = "python bundle/remote.py bundle --fetch-uma; rc=$?; rm -f .hf_token; exit $rc"
DOCTOR = "python bundle/remote.py bundle --doctor"


def hf_token():
    """手元の HF トークン（環境変数 HF_TOKEN か `hf auth login` で保存したもの）。"""
    tok = os.environ.get("HF_TOKEN")
    p = Path(os.environ.get("HF_HOME", Path.home() / ".cache/huggingface")) / "token"
    if not tok and p.exists():
        tok = p.read_text().strip()
    if not tok:
        raise SystemExit("HF のトークンがありません。`hf auth login` を実行してください")
    return tok


class Run:
    def __init__(self, ini, yes, keep):
        self.S = local.load_settings(ini)
        self.yes, self.keep = yes, keep
        self.wd = self.S["workdir"]
        self.work = self.wd / "work"
        (self.wd / "log").mkdir(parents=True, exist_ok=True)
        self.logf = self.wd / "log" / "tsrelay.log"

    def log(self, msg):
        line = f"{datetime.now():%Y-%m-%d %H:%M:%S} {msg}"
        print(line, flush=True)
        with open(self.logf, "a") as f:
            f.write(line + "\n")

    def stage(self, name):
        d = self.work / name
        return d, (d / "done").exists()

    def keep_settings(self):
        """settings.ini を workdir に保存する。前回と違えば差分を出して止まる。"""
        text = self.S["ini"].read_text()
        saved = self.wd / "settings.ini"
        if saved.exists() and saved.resolve() != self.S["ini"] and saved.read_text() != text:
            diff = "".join(difflib.unified_diff(saved.read_text().splitlines(True), text.splitlines(True),
                                                "前回", "今回"))
            raise SystemExit(f"settings.ini が前回と違います。新しい workdir で回してください:\n{diff}")
        if not saved.exists():
            saved.write_text(text)

    def main(self):
        S = self.S
        self.keep_settings()
        sub = local.read_xyz(S["substrate"])[0][0]
        cat = local.read_xyz(S["catalyst"])[0][0]
        local.validate(S, sub, cat)
        core, metals = local.detect_core(S, sub, cat)

        d, done = self.stage("00_check")
        if not done:
            local.confirm(S, sub, cat, core, ask=not self.yes)
            d.mkdir(parents=True, exist_ok=True)
            (d / "core.json").write_text(json.dumps(dict(core=core, metals=metals)))
            (d / "done").touch()

        d1, done = self.stage("01_dock")
        if not done:
            shutil.rmtree(d1, ignore_errors=True)
            d1.mkdir(parents=True)
            self.log("xtb dock を実行します")
            t0 = time.time()
            local.run_dock(S, d1, core)
            self.log(f"xtb dock が終わりました（{time.time() - t0:.0f} 秒）")
            local.collect_candidates(S, d1, sub, cat, core, metals, self.log)
            (d1 / "done").touch()
        cands = json.loads((d1 / "candidates.json").read_text())

        d2, done = self.stage("02_pick")
        if not done:
            ids = pick.serve(d1, cands, core, S["directed"])
            if not ids:
                self.log("何も選ばれなかったので終了します（再実行で選択画面をもう一度出します）")
                return
            d2.mkdir(parents=True, exist_ok=True)
            (d2 / "selected.json").write_text(json.dumps(ids))
            (d2 / "done").touch()
            self.log(f"選択: {', '.join(ids)}")
        selected = json.loads((d2 / "selected.json").read_text())

        d3, done = self.stage("03_remote")
        out = self.wd / "out"
        if not done:
            todo = [i for i in selected if not (out / i / "result.json").exists()]
            if todo:
                d3.mkdir(parents=True, exist_ok=True)
                self.remote(d1, d3, out, todo, metals, len(sub))
            if all((out / i / "result.json").exists() for i in selected):
                (d3 / "done").touch()
        self.summary(selected, out)

    # ---------- リモート段 ----------

    def remote(self, d1, d3, out, todo, metals, n_sub):
        S, s = self.S, self.S["session"]
        for attempt in range(3):
            try:
                colab.ensure_session(s, "T4", self.log)
                break
            except colab.ColabError as e:
                if e.auth:
                    raise SystemExit("colab の認証が通りません。`colab --auth oauth2 sessions` を実行して確かめてください")
                self.log(f"T4 のランタイムを作れません（{attempt + 1}/3）: {str(e)[:300]}")
                if attempt == 2:
                    raise SystemExit("T4 が割り当てられないため終了します。同じコマンドで続きから再開できます")
                time.sleep(1800)
        try:
            self.prepare(d1, d3, todo, metals, n_sub)
            self.compute(d3, out)
        except colab.ColabError as e:
            if e.lost:
                self.log("ランタイムが消えました。回収済みの結果は残っています。同じコマンドで続きから再開できます")
                sys.exit(2)
            if not self.keep:
                colab.stop(s)
            raise
        if self.keep:
            self.log(f"{s}: ランタイムは残しています（止めるときは colab stop -s {s}）")
        else:
            colab.stop(s)
            self.log(f"{s}: ランタイムを止めました")

    def prepare(self, d1, d3, todo, metals, n_sub):
        S, s = self.S, self.S["session"]
        job = dict(charge=S["charge"], mult=S["mult"], bonds=S["bonds"], nstep=S["nstep"],
                   fmax_opt=S["fmax_opt"], fmax_ts=S["fmax_ts"], metals=metals, n_sub=n_sub,
                   checkpoint=f"{colab.RDIR}/{S['model']}.pt", uma_file=f"checkpoints/{S['model']}.pt")
        (d3 / "job.json").write_text(json.dumps(job))
        bundle = d3 / "bundle.tar.gz"
        with tarfile.open(bundle, "w:gz") as tar:
            tar.add(HERE / "remote.py", arcname="remote.py")
            tar.add(d3 / "job.json", arcname="job.json")
            for i in todo:
                tar.add(d1 / f"{i}.xyz", arcname=f"cands/{i}.xyz")
        colab.upload(s, bundle, f"{colab.RDIR}/bundle.tar.gz")
        st = colab.py(s, f"""
import shutil, tarfile
shutil.rmtree("bundle", ignore_errors=True)
with tarfile.open("bundle.tar.gz") as t:
    t.extractall("bundle", filter="data")
print("TSRELAY_JSON " + json.dumps({{"uma": os.path.exists({job['checkpoint']!r})}}))
""", retries=3)
        self.log(f"{len(todo)} 構造を送りました。環境を準備します")
        self.detached(s, "setup", SETUP, "Colab の環境の準備")
        if not st["uma"]:
            # 重みは手元に置かず、ランタイムが HF から直接取る（取得が済んでから置くので、あれば完全なもの）
            self.log(f"UMA の重み {S['model']} を Hugging Face から取得します")
            colab.upload_secret(s, hf_token(), f"{colab.RDIR}/.hf_token")
            try:
                self.detached(s, "uma", FETCH, "UMA の重みの取得")
            except colab.ColabError:
                # 起動できなかったときもトークンのファイルを残さない
                colab.py(s, 'os.path.exists(".hf_token") and os.remove(".hf_token")\nprint("TSRELAY_JSON {}")',
                         retries=3)
                raise
        self.detached(s, "doctor", DOCTOR, "doctor")
        self.log("doctor が通りました")

    def detached(self, s, tag, cmd, what, timeout=1800):
        """cmd を切り離して起動し、終わるまで待つ。rc≠0 なら止まる。"""
        colab.start_detached(s, tag, cmd)
        t0 = time.time()
        while (info := colab.peek(s, tag))["rc"] is None:
            if time.time() - t0 > timeout:
                raise SystemExit(f"{what}が {timeout} 秒で終わりませんでした")
            time.sleep(15)
        if info["rc"] != 0:
            raise SystemExit(f"{what}が失敗しました（rc={info['rc']}）:\n{info['tail']}")

    def compute(self, d3, out):
        s = self.S["session"]
        colab.start_detached(s, "run", "python bundle/remote.py bundle")
        self.log("本計算を始めました")
        extra = 'from pathlib import Path\ninfo = {"done": sorted(p.parent.name for p in Path("bundle/out").glob("*/done"))}'
        fetched, last = set(), ""
        while True:
            info = colab.peek(s, "run", extra)
            new = [i for i in info["done"] if i not in fetched]
            if new:
                self.fetch(s, d3, out, new)
                fetched |= set(new)
            tail = info["tail"].strip().splitlines()[-1:] or [""]
            if tail[0] != last:
                self.log(f"[remote] {tail[0]}")
                last = tail[0]
            if info["rc"] is not None:
                if info["rc"] != 0:
                    self.log(f"remote.py が rc={info['rc']} で終わりました:\n{info['tail']}")
                return
            time.sleep(POLL_S)

    def fetch(self, s, d3, out, ids):
        """済んだ構造をその都度回収する（ランタイムが消えても済んだ分は残る）。"""
        colab.py(s, f"""
import tarfile
with tarfile.open("fetch.tar.gz", "w:gz") as t:
    for i in {ids!r}:
        t.add(f"bundle/out/{{i}}", arcname=i)
print("TSRELAY_JSON {{}}")
""", timeout=300)
        local_tar = d3 / "fetch.tar.gz"
        colab.download(s, f"{colab.RDIR}/fetch.tar.gz", local_tar)
        out.mkdir(parents=True, exist_ok=True)
        for i in ids:
            shutil.rmtree(out / i, ignore_errors=True)
        with tarfile.open(local_tar) as tar:
            tar.extractall(out, filter="data")
        for i in ids:
            r = json.loads((out / i / "result.json").read_text())
            self.log(f"回収: {i} {r['status']} {r.get('reason', '')}")

    # ---------- 出力 ----------

    def summary(self, selected, out):
        rows = []
        for i in selected:
            p = out / i / "result.json"
            if not p.exists():
                rows.append(dict(id=i, status="pending"))
                continue
            r = json.loads(p.read_text())
            status = r["status"] if r["status"] == "ok" else f"failed: {r.get('reason', '')}"
            if r.get("warnings"):
                status += "（" + " / ".join(r["warnings"]) + "）"
            rows.append(dict(id=i, status=status, E_sm_eV=r.get("E_sm"), E_ts_eV=r.get("E_ts"), E_p_eV=r.get("E_p"),
                             dE_kcal=None if r.get("dE_kcal") is None else round(r["dE_kcal"], 2),
                             imag_cm=";".join(f"{f:.1f}" for f in r.get("imag_freqs", []))))
        rows.sort(key=lambda r: (not r["status"].startswith("ok"),
                                 r["dE_kcal"] if r.get("dE_kcal") is not None else 1e9))
        cols = ["id", "status", "E_sm_eV", "E_ts_eV", "E_p_eV", "dE_kcal", "imag_cm"]
        out.mkdir(parents=True, exist_ok=True)
        with open(out / "summary.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        for r in rows:
            self.log(f"{r['id']}: ΔE‡ = {r.get('dE_kcal')} kcal/mol, 虚数 {r.get('imag_cm', '')} cm⁻¹, {r['status']}")


def main():
    ap = argparse.ArgumentParser(description="xtb dock から UMA の TS 振動解析までを自走させる")
    ap.add_argument("settings")
    ap.add_argument("--yes", action="store_true", help="起動時の y/n 確認を省く")
    ap.add_argument("--keep", action="store_true", help="終わってもランタイムを止めない（調べるとき用）")
    a = ap.parse_args()
    Run(a.settings, a.yes, a.keep).main()


if __name__ == "__main__":
    main()
