"""google-colab-cli の薄いラッパー。目印行の読み取りと障害の判定もここで行う。

ランタイム上の作業場所は RDIR。返したい値はランタイム側で最後に「TSRELAY_JSON {...}」の 1 行で print させ、
手元はその行が来た時点で colab のプロセスを止めて結果とする（colab は終わっても戻らないことがあるため）。
"""
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
import threading
import time
from pathlib import Path

MARK = "TSRELAY_JSON "
RDIR = "/content/tsrelay"
AUTH_PAT = re.compile(r"auth|credential|login|unauthori[sz]ed|\b401\b|\b403\b|token", re.I)
LOST_PAT = re.compile(r"appears to be lost", re.I)
TOKEN_PAT = re.compile(r"(token=)[^&\s)']+")  # エラー文に出るランタイムのプロキシトークンはログに残さない


class ColabError(RuntimeError):
    def __init__(self, msg, auth=False, lost=False):
        super().__init__(msg)
        self.auth, self.lost = auth, lost


def _colab():
    exe = shutil.which("colab") or str(Path.home() / ".local/bin/colab")
    if not Path(exe).exists():
        raise SystemExit("colab が見つかりません。`uv tool install google-colab-cli` で入れてください")
    return exe


def run(args, stdin=None, timeout=120, check=True, done=None):
    """colab を実行して出力を返す。出力に done で始まる行が来たら、その時点で打ち切る。"""
    cmd = [_colab(), "--auth", "oauth2", *args]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    lines, seen = [], threading.Event()

    def read():
        for line in p.stdout:
            lines.append(line)
            if done and line.startswith(done):  # 送ったコードの中の文字列には反応しない
                seen.set()

    reader = threading.Thread(target=read, daemon=True)
    reader.start()
    if stdin:
        p.stdin.write(stdin)
    p.stdin.close()
    t0 = time.time()
    while p.poll() is None and not seen.is_set():
        if time.time() - t0 > timeout:
            p.kill()
            raise ColabError(f"timeout: {' '.join(args[:2])}")
        time.sleep(0.5)
    if seen.is_set() and p.poll() is None:
        p.kill()
    reader.join(5)
    out = TOKEN_PAT.sub(r"\1REDACTED", "".join(lines))
    if check and not seen.is_set() and p.returncode != 0:
        lost = bool(LOST_PAT.search(out))  # 消失は "(404/401)" を含んでも認証切れではない
        raise ColabError(f"colab {' '.join(args[:2])} → exit {p.returncode}\n{out.strip()[-2000:]}",
                         auth=bool(AUTH_PAT.search(out)) and not lost, lost=lost)
    return out


def auth_ok():
    try:
        run(["sessions"], timeout=120)
        return True
    except ColabError:
        return False


def exists(session):
    return f"[{session}]" in run(["sessions"], timeout=120, check=False)


def ensure_session(session, gpu="T4", log=print):
    """ランタイムを用意する。colab new が戻らなくても、一覧に出ていれば作成済みとみなす。"""
    if not exists(session):
        try:
            run(["new", "-s", session, "--gpu", gpu], timeout=600)
        except ColabError:
            if not exists(session):
                raise
            log(f"{session}: colab new が戻りませんでしたが、ランタイムは作成済みです")
    for attempt in range(3):  # 作った直後は最初の exec が戻らないことがある
        try:
            py(session, 'print("TSRELAY_JSON {}")', timeout=120)
            return
        except ColabError as e:
            if e.lost or attempt == 2:
                raise
            log(f"{session}: 応答がないので待ちます（{attempt + 1}/3）")
            time.sleep(30)


def _retry(fn, retries):
    """exec や転送はときどき戻らなくなるので、何度やっても同じ結果になる操作だけ再試行する。"""
    for attempt in range(retries):
        try:
            return fn()
        except ColabError as e:
            if e.lost or e.auth or attempt == retries - 1:
                raise
            time.sleep(10)


def py(session, body, timeout=60, retries=1):
    """RDIR でコードを実行し、MARK 付きで print された JSON を返す。retries は読むだけの処理にだけ使う。"""
    code = f"import json, os, sys, subprocess\nos.makedirs({RDIR!r}, exist_ok=True)\nos.chdir({RDIR!r})\n{body}\n"

    def once():
        out = run(["exec", "-s", session, "--timeout", str(timeout)], stdin=code, timeout=timeout + 60, done=MARK)
        for line in reversed(out.splitlines()):
            if line.startswith(MARK):
                return json.loads(line[len(MARK):])
        raise ColabError(f"ランタイムからの応答を読めません（{session}）:\n{out.strip()[-2000:]}",
                         lost=bool(LOST_PAT.search(out)))
    return _retry(once, retries)


def upload(session, local, remote, timeout=600):
    _retry(lambda: run(["upload", "-s", session, str(local), remote], timeout=timeout), 3)


def upload_secret(session, text, remote):
    """秘密の文字列を小さなファイルとして送る。exec のコードに入れると履歴に残るため。"""
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "secret"
        p.touch(mode=0o600)
        p.write_text(text)
        upload(session, p, remote)


def download(session, remote, local, timeout=1800):
    _retry(lambda: run(["download", "-s", session, remote, str(local)], timeout=timeout), 3)


def start_detached(session, tag, cmd):
    """cmd を切り離して起動する。出力は tag.log、終了コードは tag.rc に書く。

    exec が戻らなくても起動はできていることがあるので、ログの有無で確かめてから起動し直す。
    """
    body = f"""
for f in ({tag + '.rc'!r}, {tag + '.log'!r}):
    if os.path.exists(f):
        os.remove(f)
subprocess.Popen(["bash", "-c", "(" + {cmd!r} + ") > {tag}.log 2>&1; echo $? > {tag}.rc"], start_new_session=True,
                 stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
print("TSRELAY_JSON {{}}")
"""
    for attempt in range(3):
        try:
            py(session, body)
            return
        except ColabError as e:
            if e.lost or e.auth or attempt == 2:
                raise
            if peek(session, tag)["started"]:
                return


def peek(session, tag, extra="info = {}"):
    """tag の終了コード（未終了なら None）とログの末尾を返す。extra で info に項目を足せる。"""
    return py(session, f"""
{extra}
info["started"] = os.path.exists({tag + '.log'!r})
info["rc"] = int(open({tag + '.rc'!r}).read()) if os.path.exists({tag + '.rc'!r}) else None
info["tail"] = open({tag + '.log'!r}).read()[-1500:] if os.path.exists({tag + '.log'!r}) else ""
print("TSRELAY_JSON " + json.dumps(info))
""", retries=3)


def stop(session):
    """ランタイムを止める。既に消えていても失敗扱いにしない。"""
    run(["stop", "-s", session], timeout=300, check=False)


def safe_extract(archive, dest):
    with tarfile.open(archive) as tar:
        tar.extractall(dest, filter="data")
