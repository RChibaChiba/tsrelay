"""ローカル段: settings.ini の読込と検証、コア検出、xtb dock、候補の後処理。

ini に書く番号は 1 始まり。0 始まりへの変換は load_settings だけで行う。
錯体の番号は「基質が先、触媒が後」（触媒の k 番目 = 基質の原子数 + k）。
"""
import configparser
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import covalent_radii

HARTREE_KCAL = 627.509474
# 3〜12 族、ランタノイド、アクチノイド、Al・Ga・In・Sn・Pb・Bi
METALS = set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81)) | set(range(89, 113)) | {13, 31, 49, 50, 82, 83}
DIRECTED_SCALE = 1.5  # directed 原子と金属の距離の上限（共有結合半径の和に掛ける）
CORE_TOL = 0.05       # 重ね合わせ後のコアのずれの上限（Å）
DUP_RMSD = 0.25       # 重複とみなす基質の重原子の RMSD（Å）


# ---------- xyz とジオメトリの小道具 ----------

def read_xyz(path):
    """複数構造の xyz を読み、(Atoms, コメント行) のリストを返す。"""
    lines = Path(path).read_text().splitlines()
    out, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        rows = [l.split() for l in lines[i + 2:i + 2 + n]]
        atoms = Atoms([r[0] for r in rows], positions=[[float(x) for x in r[1:4]] for r in rows])
        out.append((atoms, lines[i + 1]))
        i += n + 2
    return out


def write_xyz(path, atoms, comment=""):
    rows = [f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}" for s, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.positions)]
    Path(path).write_text(f"{len(atoms)}\n{comment}\n" + "\n".join(rows) + "\n")


def fit(P, Q):
    """P @ R + t が Q に最も重なる回転 R と並進 t（Kabsch）。"""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    return R, qc - pc @ R


def bonded(Z, X, i, j, scale):
    return np.linalg.norm(X[i] - X[j]) <= (covalent_radii[Z[i]] + covalent_radii[Z[j]]) * scale


def bond_set(atoms, idx, scale, skip=()):
    """idx の中の結合の組。skip に含まれる原子が関わる結合は除く。"""
    Z, X = atoms.numbers, atoms.positions
    idx = [i for i in idx if i not in skip]
    return {(i, j) for a, i in enumerate(idx) for j in idx[a + 1:] if bonded(Z, X, i, j, scale)}


# ---------- 設定 ----------

def _ints(text):
    """"5,6" や "1-3,7" を整数のリストにする。"""
    out = []
    for part in (text or "").replace(" ", "").split(","):
        if "-" in part:
            a, b = map(int, part.split("-"))
            out += range(a, b + 1)
        elif part:
            out.append(int(part))
    return out


def load_settings(path):
    path = Path(path).resolve()
    c = configparser.ConfigParser()
    c.read_string(path.read_text())
    base = path.parent

    def get(sec, key, fallback=None):
        v = c.get(sec, key, fallback=fallback)
        if v is None:
            raise SystemExit(f"settings.ini に [{sec}] {key} がありません")
        return v

    bonds = []
    for item in get("scan", "bonds").split(","):
        pair, dist = item.split(":")
        i, j = map(int, pair.split("-"))
        bonds.append((i - 1, j - 1, float(dist)))
    return dict(
        ini=path,
        substrate=base / get("input", "substrate"),
        catalyst=base / get("input", "catalyst"),
        workdir=base / get("input", "workdir"),
        charge=int(get("state", "total_charge")),
        mult=int(get("state", "total_multiplicity")),
        sub_charge=int(get("state", "substrate_charge")),
        sub_mult=int(get("state", "substrate_multiplicity")),
        cat_charge=int(get("state", "catalyst_charge")),
        cat_mult=int(get("state", "catalyst_multiplicity")),
        directed=[i - 1 for i in _ints(get("dock", "directed_atoms"))],
        xtb_args=get("dock", "xtb_args", "").split(),
        cutoff_scale=float(get("core", "cutoff_scale", "1.2")),
        add_atoms=[i - 1 for i in _ints(get("core", "add_atoms", ""))],
        drop_atoms=[i - 1 for i in _ints(get("core", "drop_atoms", ""))],
        bonds=bonds,
        nstep=int(get("scan", "nstep", "10")),
        checkpoint=Path(os.path.expanduser(get("uma", "checkpoint"))),
        fmax_opt=float(get("ts", "fmax_opt", "0.02")),
        fmax_ts=float(get("ts", "fmax_ts", "0.005")),
        cores=int(get("run", "cores", str(os.cpu_count()))),
        session=get("run", "session", "tsrelay"),
    )


def validate(S, sub, cat):
    """計算を始める前の検査。不成立なら理由を並べて止まる。"""
    errs = []
    if S["sub_charge"] + S["cat_charge"] != S["charge"]:
        errs.append(f"電荷の和が合いません: 基質 {S['sub_charge']} + 触媒 {S['cat_charge']} ≠ 全体 {S['charge']}")
    zsum = int(sub.numbers.sum() + cat.numbers.sum())
    for name, z, q, m in [("基質", sub.numbers.sum(), S["sub_charge"], S["sub_mult"]),
                          ("触媒", cat.numbers.sum(), S["cat_charge"], S["cat_mult"]),
                          ("全体", zsum, S["charge"], S["mult"])]:
        ne = int(z) - q
        if m < 1 or ne < 0 or (ne - (m - 1)) % 2:
            errs.append(f"{name}: 電子数 {ne} と多重度 {m} の偶奇が合いません")
    n_sub, n = len(sub), len(sub) + len(cat)
    if not S["directed"]:
        errs.append("directed_atoms が空です")
    for i in S["directed"]:
        if not 0 <= i < n_sub:
            errs.append(f"directed_atoms の {i + 1} が基質の原子数 {n_sub} を超えています")
    for i, j, d in S["bonds"]:
        if not (0 <= i < n and 0 <= j < n) or i == j or d <= 0:
            errs.append(f"scan.bonds の {i + 1}-{j + 1}:{d} が不正です（錯体の原子数 {n}）")
    for i in S["add_atoms"] + S["drop_atoms"]:
        if not 0 <= i < len(cat):
            errs.append(f"core の番号 {i + 1} が触媒の原子数 {len(cat)} を超えています")
    if errs:
        raise SystemExit("設定エラー:\n  " + "\n  ".join(errs))


def write_chrg_uhf(S, d):
    """xtb dock 用の .CHRG と .UHF（上から全体、input1、input2）。多重度 → 不対電子数の変換はここだけ。"""
    Path(d, ".CHRG").write_text(f"{S['charge']}\n{S['sub_charge']}\n{S['cat_charge']}\n")
    Path(d, ".UHF").write_text(f"{S['mult'] - 1}\n{S['sub_mult'] - 1}\n{S['cat_mult'] - 1}\n")


# ---------- コア ----------

def detect_core(S, sub, cat):
    """金属 + 配位原子を錯体の番号で返す。金属が複数あればすべてを対象にする。"""
    Z, X = cat.numbers, cat.positions
    metals = [k for k, z in enumerate(Z) if z in METALS]
    if not metals:
        raise SystemExit("触媒に金属が見つかりません")
    core = set(metals)
    for m in metals:
        core |= {k for k in range(len(cat)) if k != m and bonded(Z, X, m, k, S["cutoff_scale"])}
    core = (core | set(S["add_atoms"])) - set(S["drop_atoms"])
    n = len(sub)
    return sorted(n + k for k in core), [n + m for m in metals]


def confirm(S, sub, cat, core, ask=True):
    """コア原子とスキャンするペアを表示し、y/n で確かめる。"""
    cx = sub + cat
    sym = cx.get_chemical_symbols()
    n_sub = len(sub)
    print("コア（錯体の番号）: " + ", ".join(f"{sym[i]}{i + 1}" for i in core))
    print("directed（基質の番号）: " + ", ".join(f"{sym[i]}{i + 1}" for i in S["directed"]))
    for i, j, d in S["bonds"]:
        same = (i < n_sub) == (j < n_sub)
        now = f"現在 {cx.get_distance(i, j):.2f} Å" if same else "現在の距離は dock 後に決まる"
        print(f"scan: {sym[i]}{i + 1}-{sym[j]}{j + 1} → {d:.2f} Å（{now}）")
    if ask and input("この内容で進めますか [y/N]: ").strip().lower() != "y":
        raise SystemExit("中止しました")


# ---------- xtb dock ----------

def _xtb_env(S):
    xtb = os.environ.get("XTB") or shutil.which("xtb")
    if not xtb:
        raise SystemExit("xtb が見つかりません")
    env = dict(os.environ, OMP_NUM_THREADS=str(S["cores"]), OMP_STACKSIZE="1G")
    share = Path(os.path.realpath(xtb)).parent.parent / "share" / "xtb"
    if "XTBPATH" not in env and share.is_dir():
        env["XTBPATH"] = str(share)
    return xtb, env


def run_dock(S, d, core):
    d = Path(d)
    shutil.copy(S["substrate"], d / "substrate.xyz")
    shutil.copy(S["catalyst"], d / "catalyst.xyz")
    write_chrg_uhf(S, d)
    one = lambda idx: ",".join(str(i + 1) for i in idx)
    (d / "xcontrol.inp").write_text(
        f"$fix\n   atoms: {one(core)}\n$end\n$directed\n   atoms: {one(S['directed'])}\n$end\n")
    xtb, env = _xtb_env(S)
    cmd = [xtb, "dock", "substrate.xyz", "catalyst.xyz", "--input", "xcontrol.inp", *S["xtb_args"]]
    with open(d / "dock.log", "w") as log:
        rc = subprocess.run(["bash", "-c", 'ulimit -s unlimited; exec "$@"', "_", *cmd],
                            cwd=d, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    out = d / "final_structures.xyz"
    if rc != 0 or not out.exists() or out.stat().st_size == 0:
        raise SystemExit(f"xtb dock が失敗しました（rc={rc}）。{d / 'dock.log'} を見てください")


def collect_candidates(S, d, sub, cat, core, metals, log):
    """dock の構造を検査し、残ったものを cand_XXX.xyz と candidates.json に書く。"""
    d = Path(d)
    n_sub = len(sub)
    ref = sub + cat
    scale = S["cutoff_scale"]
    Zc = ref.numbers
    sub_idx, cat_idx = range(n_sub), range(n_sub, len(ref))
    ref_bonds = bond_set(ref, sub_idx, scale, metals) | bond_set(ref, cat_idx, scale, metals)
    heavy = [i for i in sub_idx if Zc[i] > 1] or list(sub_idx)
    kept = []  # (Atoms をコアで入力に重ねたもの, 情報)
    for k, (atoms, comment) in enumerate(read_xyz(d / "final_structures.xyz")):
        cid = f"cand_{k:03d}"
        e = float(comment.split()[0])
        R, t = fit(atoms.positions[core], ref.positions[core])
        X = atoms.positions @ R + t
        dev = np.linalg.norm(X[core] - ref.positions[core], axis=1).max()
        atoms.positions = X
        dmin = [min(atoms.get_distance(i, m) for m in metals) for i in S["directed"]]
        far = [i for i in S["directed"]
               if min(atoms.get_distance(i, m) / (covalent_radii[Zc[i]] + covalent_radii[Zc[m]]) for m in metals)
               > DIRECTED_SCALE]
        bonds = bond_set(atoms, sub_idx, scale, metals) | bond_set(atoms, cat_idx, scale, metals)
        dup = next((c["id"] for a, c in kept
                    if np.sqrt(((a.positions[heavy] - X[heavy]) ** 2).sum(1).mean()) < DUP_RMSD), None)
        if dev > CORE_TOL:
            reason = f"コアが {dev:.3f} Å ずれた"
        elif bonds != ref_bonds:
            reason = f"分子内の結合が変わった（{len(bonds ^ ref_bonds)} 本）"
        elif far:
            reason = "directed 原子が金属を向いていない（" + ",".join(str(i + 1) for i in far) + "）"
        elif dup:
            reason = f"{dup} と重複"
        else:
            kept.append((atoms, dict(id=cid, e_h=e, directed_metal=[round(x, 3) for x in dmin])))
            continue
        log(f"{cid}: 捨てた（{reason}）")
    if not kept:
        raise SystemExit("dock の候補がすべて検査で落ちました。ログを見てください")
    e0 = min(c["e_h"] for _, c in kept)
    for atoms, c in kept:
        c["rel_kcal"] = round((c["e_h"] - e0) * HARTREE_KCAL, 2)
        write_xyz(d / f"{c['id']}.xyz", atoms, f"{c['id']} E={c['e_h']:.8f}")
    cands = [c for _, c in kept]
    (d / "candidates.json").write_text(json.dumps(cands, indent=1))
    log(f"dock 候補: {len(cands)} 個を残した")
    return cands
