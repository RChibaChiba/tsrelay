"""リモート段（Colab T4）: 選んだ構造ごとに UMA で p 錯体を整え、逆スキャンで sm 錯体を作り、DMF と Sella で TS を求める。

    python remote.py BUNDLE_DIR            # 本計算
    python remote.py BUNDLE_DIR --doctor   # 環境の検査
    python remote.py BUNDLE_DIR --fetch-uma  # UMA の重みを HF から取る

BUNDLE_DIR には job.json と cands/<ID>.xyz を置く。結果は out/<ID>/ に書き、1 構造済むごとに done を置く。
番号はすべて 0 始まり（変換は手元の local.load_settings で済んでいる）。
"""
import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from ase.data import covalent_radii
from ase.io import read, write

# 版は固定する。変えたら Colab で doctor を通してから
PINS = {"fairchem-core": "2.25.0", "pydmf": "1.2.2", "sella": "2.6.0", "cyipopt": "1.7.0"}
EV_KCAL = 23.060548
BOND_BACK = 1.3     # 伸ばしたペアがこれ × 共有結合半径の和より長ければ「結合は戻っていない」
SM_SHIFT = 1.0      # 基質–金属の最短距離が p からこれ以上変わったら警告（Å）
SMALL_IMAG = 50.0   # これ未満の虚数振動は数えない（cm⁻¹）
NOISE_IMAG = 10.0   # これ未満は並進・回転の数値ノイズとして注記もしない（cm⁻¹）
MODE_FRAC = 0.2     # 虚数モードのうち scan 原子が占める割合がこれ未満なら警告
NMOVE = 10          # DMF の可動点の数
FD_STEP = 0.005     # 有限差分ヘッセの変位（Å）

UNIT = None
JOB = None


def fit(P, Q):
    """P @ R + t が Q に最も重なる回転 R と並進 t（Kabsch）。"""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    return R, qc - pc @ R


# ---------- UMA ----------

def load_unit(checkpoint):
    import torch
    from fairchem.core import pretrained_mlip
    device = "cuda" if torch.cuda.is_available() else "cpu"
    return pretrained_mlip.load_predict_unit(str(checkpoint), device=device)


def setup(atoms):
    from fairchem.core import FAIRChemCalculator
    atoms.info.update(charge=JOB["charge"], spin=JOB["mult"])
    atoms.calc = FAIRChemCalculator(UNIT, task_name="omol")
    return atoms


def bare(atoms):
    """計算器と拘束を外した複製。"""
    a = atoms.copy()
    a.set_constraint()
    a.calc = None
    return a


def energy(atoms):
    return float(atoms.get_potential_energy())


# ---------- 最適化と解析 ----------

def minimize(atoms, fmax, steps=2000):
    """拘束なしの最小化。Sella（内部座標）で回し、内部座標が作れないときだけ LBFGS に落とす。"""
    from sella import Sella
    try:
        return Sella(atoms, order=0, internal=True, logfile=None).run(fmax=fmax, steps=steps)
    except Exception:
        from ase.optimize import LBFGS
        return LBFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)


def reverse_scan(atoms, bonds, nstep, fmax):
    """各ペアを現在の距離から目標距離まで nstep 刻みで動かす。各刻みで全ペアを固定して最適化する。"""
    from ase.constraints import FixBondLengths
    from ase.optimize import LBFGS
    d0 = [atoms.get_distance(i, j) for i, j, _ in bonds]
    pairs = [(i, j) for i, j, _ in bonds]
    frames = []
    for k in range(1, nstep + 1):
        targets = [a + (t - a) * k / nstep for a, (_, _, t) in zip(d0, bonds)]
        atoms.set_constraint()
        for _ in range(3):  # ペアが原子を共有していても目標値に寄せる
            for (i, j), d in zip(pairs, targets):
                atoms.set_distance(i, j, d, fix=0.5)
        # 既定の tolerance=1e-13 では UMA の力で収束判定が通らず例外になる
        atoms.set_constraint(FixBondLengths(pairs, bondlengths=targets, tolerance=1e-6, iterations=1000))
        LBFGS(atoms, logfile=None).run(fmax=fmax, steps=500)
        frames.append(bare(atoms))
    atoms.set_constraint()
    return frames


def min_sub_metal(atoms):
    X, n = atoms.positions, JOB["n_sub"]
    return min(np.linalg.norm(X[i] - X[m]) for i in range(n) for m in JOB["metals"])


def check_sm(sm, p):
    """伸ばしたペアが戻っていないか（失敗）と、基質–金属の距離が大きく変わっていないか（警告）。

    縮めたペア（目標距離が p より短いもの）は結合を作る側なので、戻りの判定には使わない。
    """
    Z = sm.numbers
    for i, j, target in JOB["bonds"]:
        if target <= p.get_distance(i, j):
            continue
        d = sm.get_distance(i, j)
        lim = BOND_BACK * (covalent_radii[Z[i]] + covalent_radii[Z[j]])
        if d <= lim:
            return f"sm で {i + 1}-{j + 1} が {d:.2f} Å に戻った（{lim:.2f} Å 以下）", None
    shift = min_sub_metal(sm) - min_sub_metal(p)
    warn = f"基質–金属の最短距離が p から {shift:+.2f} Å 変わった" if abs(shift) >= SM_SHIFT else None
    return None, warn


def run_dmf(sm, p, out):
    """sm を p に重ねてから DMF で経路を求め、最高エネルギー点の推定を返す。"""
    from dmf import DirectMaxFlux, interpolate_fbenm
    sm = bare(sm)
    R, t = fit(sm.positions, p.positions)
    sm.positions = sm.positions @ R + t
    ref = [sm, bare(p)]
    fb = interpolate_fbenm(ref, correlated=True)
    mx = DirectMaxFlux(ref, coefs=fb.coefs.copy(), nmove=NMOVE, update_teval=True)
    for im in mx.images:
        setup(im)
    mx.add_ipopt_options({"output_file": str(out / "dmf_ipopt.out")})
    mx.solve(tol="middle")
    write(out / "path.xyz", [bare(im) for im in mx.images])
    return bare(mx.history.images_tmax[-1])


def hessian(atoms):
    """力の中心差分によるヘッセ行列（eV/Å²）。"""
    x0 = atoms.positions.copy()
    n = len(atoms)
    H = np.zeros((3 * n, 3 * n))
    for a in range(n):
        for c in range(3):
            g = []
            for s in (1, -1):
                x = x0.copy()
                x[a, c] += s * FD_STEP
                atoms.positions = x
                g.append(atoms.get_forces().ravel())
            H[3 * a + c] = -(g[0] - g[1]) / (2 * FD_STEP)
    atoms.positions = x0
    return (H + H.T) / 2


def frequencies(ts):
    """振動数（cm⁻¹、虚数は負で表す）、最も負のモード、そのモードで scan 原子が占める割合。

    割合の分母は重原子と scan 原子だけにする（H の揺れが大きさを占めて、本物の TS でも割合が低く出るため）。
    """
    from ase.vibrations import VibrationsData
    vib = VibrationsData.from_2d(ts, hessian(ts))
    en, modes = vib.get_energies_and_modes(all_atoms=True)
    cm = vib.get_frequencies()
    signed = np.where(np.abs(cm.imag) > np.abs(cm.real), -np.abs(cm.imag), cm.real)
    k = int(np.argmin(signed))
    v = (modes[k] ** 2).sum(1)
    scan = sorted({i for i, j, _ in JOB["bonds"]} | {j for i, j, _ in JOB["bonds"]})
    base = sorted(set(scan) | {i for i, z in enumerate(ts.numbers) if z > 1})
    return signed, modes[k], float(v[scan].sum() / v[base].sum())


# ---------- 1 構造の手順 ----------

def run_one(cid, src, out):
    res = dict(id=cid, status="failed", warnings=[])
    p = setup(read(src))
    minimize(p, JOB["fmax_opt"])
    res["E_p"] = energy(p)
    write(out / "p.xyz", bare(p))

    sm = setup(bare(p))
    write(out / "scan.xyz", reverse_scan(sm, JOB["bonds"], JOB["nstep"], JOB["fmax_opt"]))
    minimize(sm, JOB["fmax_opt"])
    res["E_sm"] = energy(sm)
    write(out / "sm.xyz", bare(sm))
    fail, warn = check_sm(sm, p)
    if warn:
        res["warnings"].append(warn)
    if fail:
        res["reason"] = fail
        return res

    from sella import Sella
    ts = setup(run_dmf(sm, p, out))
    write(out / "ts_guess.xyz", bare(ts))
    conv = Sella(ts, order=1, internal=True, hessian_function=hessian,
                 logfile=str(out / "sella.log")).run(fmax=JOB["fmax_ts"], steps=1000)
    res["E_ts"] = energy(ts)
    write(out / "ts.xyz", bare(ts))

    freqs, mode, frac = frequencies(ts)
    (out / "freq.txt").write_text("\n".join(f"{f:10.2f}" for f in freqs) + "\n")
    # 虚数モードを両向きに振った 2 構造（ビューアで動きを確かめる用）
    frames = [bare(ts) for _ in range(2)]
    for a, s in zip(frames, (1, -1)):
        a.positions += s * 0.5 * mode / np.abs(mode).max()
    write(out / "ts_mode.xyz", frames)
    imag = [float(f) for f in freqs if f <= -SMALL_IMAG]
    small = [float(f) for f in freqs if -SMALL_IMAG < f <= -NOISE_IMAG]
    res.update(imag_freqs=imag, scan_mode_frac=round(frac, 3),
               dE_kcal=(res["E_ts"] - res["E_sm"]) * EV_KCAL)
    if small:
        res["warnings"].append("小さな虚数振動: " + ", ".join(f"{f:.1f}" for f in small))
    if frac < MODE_FRAC:
        res["warnings"].append(f"虚数モードの scan 原子の割合が {frac:.2f}（別反応の TS の疑い）")
    if not conv:
        res["reason"] = "Sella が収束しなかった"
    elif len(imag) != 1:
        res["reason"] = f"虚数振動が {len(imag)} 個"
    else:
        res["status"] = "ok"
    return res


# ---------- doctor ----------

def doctor():
    """GPU が T4 か、版が PINS どおりか、UMA で水のエネルギーが出るかを確かめる。例外はそのまま rc≠0 になる。"""
    import torch
    from importlib.metadata import version
    from ase import Atoms
    import cyipopt, dmf, sella  # noqa: F401
    gpu = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    versions = {k: version(k) for k in PINS}
    global UNIT
    UNIT = load_unit(JOB["checkpoint"])
    w = setup(Atoms("OH2", positions=[[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]]))
    w.info.update(charge=0, spin=1)
    info = dict(gpu=gpu, versions=versions, water_eV=energy(w))
    print(json.dumps(info, indent=1))
    return 0 if gpu and "T4" in gpu and versions == PINS else 1


def fetch_uma(d):
    """UMA の重みを Hugging Face から取る。取り終えてから checkpoint に置くので、途中で切れても半端なファイルは残らない。

    トークンは bundle の隣の .hf_token で受け取り、読んだらすぐ消す（ディスクにもログにも残さない）。
    """
    from huggingface_hub import hf_hub_download
    tok = d.parent / ".hf_token"
    token = tok.read_text().strip()
    tok.unlink()
    dst = Path(JOB["checkpoint"])
    src = hf_hub_download("facebook/UMA", JOB["uma_file"], token=token, local_dir=dst.parent / "hf")
    os.replace(src, dst)
    print(f"{JOB['uma_file']} を取得しました")
    return 0


def main():
    global UNIT, JOB
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle")
    ap.add_argument("--doctor", action="store_true")
    ap.add_argument("--fetch-uma", action="store_true")
    a = ap.parse_args()
    d = Path(a.bundle).resolve()
    JOB = json.loads((d / "job.json").read_text())
    if a.fetch_uma:
        sys.exit(fetch_uma(d))
    if a.doctor:
        sys.exit(doctor())
    UNIT = load_unit(JOB["checkpoint"])
    for src in sorted((d / "cands").glob("*.xyz")):
        cid = src.stem
        out = d / "out" / cid
        if (out / "done").exists():
            continue
        out.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        try:
            res = run_one(cid, src, out)
        except Exception as e:
            res = dict(id=cid, status="failed", reason=f"{type(e).__name__}: {e}",
                       traceback=traceback.format_exc()[-3000:])
        res["seconds"] = round(time.time() - t0, 1)
        (out / "result.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
        (out / "done").touch()
        print(f"{cid}: {res['status']} {res.get('reason', '')} ({res['seconds']} s)", flush=True)


if __name__ == "__main__":
    main()
