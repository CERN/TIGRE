"""TIGRE against the ASTRA Toolbox: a wider speed suite (cone, parallel and fan beam; projectors, FDK/FBP, SIRT, CGLS).

Builds on tigre_vs_astra.py (keep both files in the same folder: this one imports its helpers). That script times
one geometry family (TIGRE's default cone beam, n^3, 360 angles) and produced the 3.3.0 note's tables; this one
varies the problem, so that the comparison does not rest on one shape. Groups (--groups):

    size       cone beam, n^3 (128..1024), detector of n x n pixels over 409.6 mm (it just sees the image),
               fp, bp, fdk; NumPy arrays on every GPU count asked for, GPU tensors on 1 GPU
    legacy     as the 3.3.0 runs: 0.8 mm detector pixels whatever n (at 1024 the detector is 2x the image's
               shadow, so most rays miss it). Shows whether the 3.3.0 1024^3 Siddon result depends on that
    angles     cone 512^3, 45..1440 angles, fp and bp
    shape      non-cubic problems: a clinical CBCT (256x512x512 image, 512x1024 detector), a fine detector
               (256^3 image, 1024^2 detector) and a coarse one (512^3 image, 256^2 detector); fp, bp, fdk
    offset     cone 512^3 with the detector shifted sideways by a quarter of its width (ASTRA cone_vec); fp, bp
    parallel   parallel beam n^3 (256..1024): TIGRE parallel mode against ASTRA parallel3d; fp, bp
    fan        fan beam (TIGRE: one-slice cone geometry; ASTRA: 2D fanflat), n^2 (512..4096), 360 and
               1440 angles; fp, bp, and TIGRE fdk against ASTRA FBP_CUDA
    iterative  TIGRE sirt / cgls against ASTRA SIRT3D_CUDA / CGLS3D_CUDA (cone 256^3, 512^3; NumPy arrays, every
               GPU count) and SIRT_CUDA / CGLS_CUDA (fan 1024^2); time of the whole call for --niter iterations
    latency    small problems (32^3..128^3, 100 angles): the fixed cost of a call, which matters in loops
    large      (not run by default) cone 1536^3 and 2048^3, fp and bp, NumPy arrays: problems split to fit the
               GPU. Needs ~100 GB of host memory at 2048^3

Same problem in both toolboxes: each beam type (cone, parallel, fan) gets its own calibration of ASTRA's angle
convention and detector flips against TIGRE's (as in tigre_vs_astra.py), and the shifted detector its own check of
the direction of the shift. Every TIGRE output is then compared with ASTRA's ("agree", "scale"): a check that the
same problem was timed, not an accuracy measure (the projection models differ, and so do backprojection weights
and the SIRT/CGLS details: TIGRE's CGLS projects 3 times per iteration, ASTRA's 2).

Timing as in tigre_vs_astra.py (rounds that alternate which toolbox goes first; per round the median of the
calls; TIGRE/ASTRA ratio per round, mean and 95% CI), except that the number of calls per round grows for fast
calls so that each measurement lasts at least --min-time seconds (at most 200 calls).

Usage (outside the TIGRE source tree, so that the installed TIGRE is imported):

    python tigre_vs_astra_suite.py --list                           # the cases that would run, nothing else
    python tigre_vs_astra_suite.py --check                          # the calibrations only, about a minute
    python tigre_vs_astra_suite.py --out suite.json                 # every default group; 1, 2 and 4 GPUs if there
    python tigre_vs_astra_suite.py --groups fan latency --out fan.json
    python tigre_vs_astra_suite.py --quick --gpus 1 --out win_suite.json   # smaller sizes, 2 rounds
    python tigre_vs_astra_suite.py --out suite.json --resume        # carry on after a crash or a stop
    python tigre_vs_astra_suite.py --report suite.json --csv suite.csv

The results file is written after every case. GPUs are 0..g-1 of the visible ones (CUDA_VISIBLE_DEVICES).

This is a benchmark, not a test: it is not collected by pytest and must not run in CI.
"""
import argparse
import json
import math
import os
import platform
import re
import statistics
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tigre_vs_astra as base  # noqa: E402  (agreement, ci, gpu_ids, count_gpus, gpu_name)

# (operation, TIGRE call, ASTRA call)
PAIRS = [
    ("fp", "TIGRE fp Siddon", "ASTRA fp"),
    ("fp", "TIGRE fp interpolated", "ASTRA fp"),
    ("bp", "TIGRE bp matched", "ASTRA bp"),
    ("bp", "TIGRE bp FDK", "ASTRA bp"),
    ("bp", "TIGRE bp", "ASTRA bp"),  # parallel beam: TIGRE has one backprojector there
    ("fdk", "TIGRE fdk", "ASTRA fdk"),  # fan beam: ASTRA's is FBP_CUDA
    ("sirt", "TIGRE sirt", "ASTRA sirt"),
    ("cgls", "TIGRE cgls", "ASTRA cgls"),
]
OPS = ["fp", "bp", "fdk", "sirt", "cgls"]
GROUPS = ["size", "legacy", "angles", "shape", "offset", "parallel", "fan", "iterative", "latency", "large", "paper"]
DEFAULT_GROUPS = [g for g in GROUPS if g not in ("large", "paper")]
DSO, DSD = 1000.0, 1536.0
FOV_DET = 409.6  # mm: TIGRE's default 512 x 0.8 mm detector, which just sees the 256 mm image


# ------------------------------------------------------------------------------------------------------ the cases

def case(group, beam, nvox, ndet, ddet, dvox, angles=360, offdet=(0.0, 0.0), ops=("fp", "bp"), multi_gpu=False,
         modes=("host", "device"), niter=None, label=None):
    nz, ny, nx = nvox
    nv, nu = ndet
    key = "{}: {} {}x{}x{}, det {}x{}, {} angles".format(group, beam, nz, ny, nx, nv, nu, angles)
    if any(offdet):
        key += ", det offset {:g},{:g} mm".format(*offdet)
    if niter:
        key += ", {} iterations".format(niter)
    return {"key": key, "group": group, "beam": beam, "nvox": list(nvox), "ndet": list(ndet),
            "ddet": [float(x) for x in ddet], "dvox": float(dvox), "angles": int(angles),
            "offdet": [float(x) for x in offdet], "ops": list(ops), "multi_gpu": multi_gpu, "modes": list(modes),
            "niter": niter,
            "label": label or ("{}^3".format(nx) if nz == ny == nx else "{}x{}x{}".format(nz, ny, nx))}


def cone_cube(group, n, angles=360, ops=("fp", "bp"), legacy=False, **kw):
    ddet = 0.8 if legacy else FOV_DET / n
    return case(group, "cone", (n, n, n), (n, n), (ddet, ddet), 256.0 / n, angles, ops=ops, **kw)


def fan_square(group, n, angles=360, ops=("fp", "bp", "fdk"), **kw):
    dvox = 256.0 / n
    kw.setdefault("label", "{}^2".format(n))
    return case(group, "fan", (1, n, n), (1, n), (dvox * DSD / DSO, FOV_DET / n), dvox, angles, ops=ops,
                multi_gpu=False, **kw)


def build_cases(groups, quick, niter, fan_niter):
    out = []
    if "size" in groups:
        for n in ([128, 256, 512] if quick else [128, 256, 512, 1024]):
            out.append(cone_cube("size", n, ops=("fp", "bp", "fdk"), multi_gpu=True))
    if "legacy" in groups:
        for n in ([256] if quick else [256, 1024]):
            out.append(cone_cube("legacy", n, ops=("fp", "bp"), legacy=True, label="{}^3, 0.8 mm".format(n)))
    if "angles" in groups:
        for a in ([90, 360, 1440] if quick else [45, 90, 180, 360, 720, 1440]):
            out.append(cone_cube("angles", 512, angles=a, label="{} angles".format(a)))
    if "shape" in groups:
        out.append(case("shape", "cone", (256, 512, 512), (512, 1024), (0.45, 0.4), 0.5,
                        ops=("fp", "bp", "fdk"), label="clinical CBCT"))
        out.append(case("shape", "cone", (256, 256, 256), (1024, 1024), (0.4, 0.4), 1.0,
                        ops=("fp", "bp", "fdk"), label="fine detector"))
        out.append(case("shape", "cone", (512, 512, 512), (256, 256), (1.6, 1.6), 0.5,
                        ops=("fp", "bp", "fdk"), label="coarse detector"))
    if "offset" in groups:
        out.append(case("offset", "cone", (512, 512, 512), (512, 512), (0.8, 0.8), 0.5, offdet=(0.0, FOV_DET / 4),
                        label="u offset 102.4 mm"))
    if "parallel" in groups:
        for n in ([256, 512] if quick else [256, 512, 1024]):
            out.append(case("parallel", "parallel", (n, n, n), (n, n), (256.0 / n, 256.0 / n), 256.0 / n,
                            multi_gpu=True))
    if "fan" in groups:
        for n in ([512, 1024, 2048] if quick else [512, 1024, 2048, 4096]):
            out.append(fan_square("fan", n))
        out.append(fan_square("fan", 1024, angles=1440, label="1024^2, 1440 angles"))
    if "iterative" in groups:
        for n in ([256] if quick else [256, 512]):
            out.append(cone_cube("iterative", n, ops=("sirt", "cgls"), multi_gpu=True, modes=("host",),
                                 niter=niter))
        out.append(fan_square("iterative", 1024, ops=("sirt", "cgls"), modes=("host",), niter=fan_niter,
                              label="fan 1024^2"))
    if "latency" in groups:
        for n in (32, 64, 128):
            out.append(cone_cube("latency", n, angles=100))
    if "paper" in groups:
        # the regime of the JPDC 2020 multi-GPU paper: N^3 image, N^2 detector, N angles (about 3x more kernel work per
        # byte moved than 360 angles at 1024^3), up to sizes that are split on the GPU; host arrays, 1/2/4 GPUs.
        # Cases needing more host memory than the machine has are skipped by the runner.
        for n in ([512, 1024] if quick else [512, 1024, 1536, 2048]):
            out.append(cone_cube("paper", n, angles=n, ops=("fp", "bp"), multi_gpu=True, modes=("host",),
                                 label="{}^3, {} angles".format(n, n)))
    if "large" in groups:
        for n in ([1536] if quick else [1536, 2048]):
            out.append(cone_cube("large", n, ops=("fp", "bp"), multi_gpu=True, modes=("host",)))
    return out


def plan(cases, gcounts, modes):
    """(case, mode, number of GPUs) to run, in order."""
    out = []
    for c in cases:
        for mode in ("host", "device"):
            if mode not in modes or mode not in c["modes"]:
                continue
            counts = gcounts if (mode == "host" and c["multi_gpu"]) else [1]
            for g in counts:
                out.append((c, mode, g))
    return out


def run_key(c, mode, g):
    return "{} | {} | {} GPU".format(c["key"], "NumPy" if mode == "host" else "tensors", g)


# ------------------------------------------------------------------------------------------------ geometry, data

def build_geometry(c):
    import tigre
    nvox = np.array(c["nvox"], dtype=np.int64)  # 64-bit: TIGRE <= 3.3.0's Atb crashed on 32-bit sizes at 1024^3
    if c["beam"] == "parallel":
        geo = tigre.geometry(mode="parallel", nVoxel=nvox)
        # TIGRE's parallel rays start at a plane DSO from the centre: the default (nVoxel[0], in mm) can be inside
        # the image (it was for the calibration's 256 mm, 64-voxel image). Put it outside, as for cone beam.
        geo.DSO, geo.DSD = DSO, DSD
    else:  # cone and fan (TIGRE's fan beam is a one-slice cone geometry)
        geo = tigre.geometry(mode="cone", default=True, nVoxel=nvox)
        geo.DSO, geo.DSD = DSO, DSD
    geo.nVoxel = nvox
    geo.dVoxel = np.array([c["dvox"]] * 3, dtype=np.float64)
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.nDetector = np.array(c["ndet"], dtype=np.int64)
    geo.dDetector = np.array(c["ddet"], dtype=np.float64)
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offOrigin = np.zeros(3)
    geo.offDetector = np.array(c["offdet"], dtype=np.float64)
    geo.rotDetector = np.zeros(3)
    geo.accuracy = 0.5
    return geo


def phantom(shape):
    """Gaussian blobs with no symmetry (as tigre_vs_astra.phantom, for any shape)."""
    nz, ny, nx = shape
    cs = [(np.arange(m, dtype=np.float64) + 0.5) / m - 0.5 for m in (nz, ny, nx)]
    img = np.zeros((nz, ny, nx), dtype=np.float32)
    blobs = [((0.10, -0.20, 0.15), 0.12, 1.0), ((-0.25, 0.10, -0.05), 0.08, 0.7),
             ((0.00, 0.25, 0.30), 0.05, 1.5), ((0.20, 0.00, -0.30), 0.10, 0.5)]
    for (z, y, x), w, amp in blobs:
        gz = np.ones(1) if nz == 1 else np.exp(-0.5 * ((cs[0] - z) / w) ** 2)
        gy, gx = (np.exp(-0.5 * ((cc - m) / w) ** 2) for cc, m in ((cs[1], y), (cs[2], x)))
        plane = np.outer(gy, gx).astype(np.float32)
        for k in range(nz):
            img[k] += np.float32(amp * gz[k]) * plane
    return img


def astra_setup(c, geo, angles, conv):
    """(volume geometry, projection geometry, 2 or 3): ASTRA's description of TIGRE's geo."""
    import astra
    nz, ny, nx = (int(v) for v in geo.nVoxel)
    sz, sy, sx = (float(v) for v in geo.sVoxel)
    nv, nu = (int(v) for v in geo.nDetector)
    dv, du = (float(v) for v in geo.dDetector)
    a = conv["sign"] * np.asarray(angles, dtype=np.float64) + conv["offset"]
    if c["beam"] == "fan":
        vol = astra.create_vol_geom(ny, nx, -sx / 2, sx / 2, -sy / 2, sy / 2)
        proj = astra.create_proj_geom("fanflat", du, nu, a, DSO, DSD - DSO)
        return vol, proj, 2
    vol = astra.create_vol_geom(ny, nx, nz, -sx / 2, sx / 2, -sy / 2, sy / 2, -sz / 2, sz / 2)
    if c["beam"] == "parallel":
        return vol, astra.create_proj_geom("parallel3d", du, dv, nv, nu, a), 3
    proj = astra.create_proj_geom("cone", du, dv, nv, nu, a, DSO, DSD - DSO)
    offv, offu = c["offdet"]
    if offv or offu:
        proj = astra.functions.geom_2vec(proj)
        vec = proj["Vectors"]  # per angle: source, detector centre, u pixel vector, v pixel vector
        su, sv = conv.get("off_u", 1), conv.get("off_v", 1)  # from calibrate_offset
        vec[:, 3:6] += (su * offu / du) * vec[:, 6:9] + (sv * offv / dv) * vec[:, 9:12]
    return vol, proj, 3


def a2t_proj(p, conv, dim):
    """ASTRA projections as TIGRE's (angles, v, u): a view."""
    if dim == 3:
        return base.astra_to_tigre(p, conv)
    q = p[:, None, :] if isinstance(p, np.ndarray) else p.unsqueeze(1)
    if conv["flip_u"]:
        q = q[:, :, ::-1] if isinstance(q, np.ndarray) else q.flip(2)
    return q


def t2a_proj(q, conv, dim):
    """TIGRE projections as ASTRA's, contiguous."""
    if dim == 3:
        return base.tigre_to_astra(q, conv)
    if isinstance(q, np.ndarray):
        q = q[:, 0, :]
        return np.ascontiguousarray(q[:, ::-1] if conv["flip_u"] else q)
    q = q[:, 0, :]
    return (q.flip(1) if conv["flip_u"] else q).contiguous()


def a2t_vol(v, dim):
    if dim == 3:
        return v
    return v[None] if isinstance(v, np.ndarray) else v.unsqueeze(0)


def t2a_vol(v, dim):
    return v if dim == 3 else v[0]


class Astra:
    """ASTRA's projector for one problem, in 2D or 3D, with the calls used here."""

    def __init__(self, vol, proj, dim):
        import astra
        import astra.experimental as ae
        self.astra, self.ae, self.vol, self.proj, self.dim = astra, ae, vol, proj, dim
        self.pid = astra.create_projector("cuda" if dim == 2 else "cuda3d", proj, vol)

    def fp(self, v, out):
        (self.ae.direct_FP2D if self.dim == 2 else self.ae.direct_FP3D)(self.pid, v, out)

    def bp(self, out, p):
        (self.ae.direct_BP2D if self.dim == 2 else self.ae.direct_BP3D)(self.pid, out, p)

    def algorithm(self, name, out, p, niter=1, options=None):
        """Run an ASTRA algorithm (FDK_CUDA, FBP_CUDA, SIRT3D_CUDA, ...) on host arrays linked, not copied."""
        astra = self.astra
        data = astra.data2d if self.dim == 2 else astra.data3d
        vid = data.link("-vol", self.vol, out)
        sid = data.link("-sino", self.proj, p)
        cfg = astra.astra_dict(name)
        cfg["ReconstructionDataId"] = vid
        cfg["ProjectionDataId"] = sid
        if self.dim == 2:
            cfg["ProjectorId"] = self.pid
        if options:
            cfg["option"] = dict(options)
        aid = astra.algorithm.create(cfg)
        try:
            astra.algorithm.run(aid, niter)
        finally:
            astra.algorithm.delete(aid)
            data.delete([vid, sid])
        return out

    def close(self):
        (self.astra.projector.delete if self.dim == 2 else self.astra.projector3d.delete)(self.pid)


# ---------------------------------------------------------------------------------------------------- calibration

def calib_case(beam, offdet=(0.0, 0.0)):
    if beam == "fan":
        return case("calibration", "fan", (1, 64, 64), (1, 96), (4 * DSD / DSO, 6.4), 4.0)
    if beam == "parallel":
        return case("calibration", "parallel", (64, 64, 64), (64, 96), (4.0, 4.0), 4.0)
    return case("calibration", "cone", (64, 64, 64), (64, 64 + 32 * bool(offdet[1])), (6.4, 6.4), 4.0,
                offdet=offdet)


def _astra_fp_host(c, geo, angles, conv, img):
    vol, proj, dim = astra_setup(c, geo, angles, conv)
    a = Astra(vol, proj, dim)
    try:
        nv, nu = (int(v) for v in geo.nDetector)
        out = np.empty((len(angles), nu) if dim == 2 else (nv, len(angles), nu), dtype=np.float32)
        a.fp(t2a_vol(img, dim), out)
    finally:
        a.close()
    return out, dim


def calibrate(beam):
    """ASTRA's angle convention and detector flips that match TIGRE, for one beam type."""
    import tigre
    c = calib_case(beam)
    geo = build_geometry(c)
    angles = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    img = phantom(c["nvox"])
    ref = tigre.Ax(img, geo, angles, "interpolated", gpuids=base.gpu_ids([0]))
    tried = []
    for sign in (1, -1):
        for k in range(4):
            conv0 = {"sign": sign, "offset": k * math.pi / 2, "flip_u": False, "flip_v": False}
            p, dim = _astra_fp_host(c, geo, angles, conv0, img)
            for flip_v in ((False,) if dim == 2 else (False, True)):
                for flip_u in (False, True):
                    conv = dict(conv0, flip_u=flip_u, flip_v=flip_v)
                    s, rel = base.agreement(ref, a2t_proj(p, conv, dim))
                    tried.append((rel, s, conv))
    tried.sort(key=lambda r: r[0])
    rel, s, conv = tried[0]
    # Conventions that give the same projections are not rivals (parallel beam: theta + 180 deg with u flipped is
    # theta): the next best is the first one that agrees measurably worse.
    others = [t for t in tried[1:] if t[0] > 1.01 * rel + 1e-9]
    tried = [tried[0]] + (others or tried[1:])
    print("{:<9} ASTRA angle = {:+d} * TIGRE angle + {:.0f} deg, flip u: {}, flip v: {}; agree {:.2e} (scale "
          "{:.4f}), next best {:.2e}".format(beam, conv["sign"], math.degrees(conv["offset"]), conv["flip_u"],
                                            conv["flip_v"], rel, s, tried[1][0]))
    if rel > 0.05 or tried[1][0] < 2 * rel:
        print("  WARNING: no convention clearly matches for {} beam: its 'agree' columns mean nothing.".format(beam))
    return {"convention": conv, "agree": rel, "scale": s, "next_best_agree": tried[1][0]}


def calibrate_offset(conv):
    """The direction in which ASTRA's cone_vec detector must move for TIGRE's offDetector (v, u)."""
    import tigre
    c = calib_case("cone", offdet=(12.8, 51.2))
    geo = build_geometry(c)
    angles = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    img = phantom(c["nvox"])
    ref = tigre.Ax(img, geo, angles, "interpolated", gpuids=base.gpu_ids([0]))
    tried = []
    for su in (1, -1):
        for sv in (1, -1):
            cv = dict(conv, off_u=su, off_v=sv)
            p, dim = _astra_fp_host(c, geo, angles, cv, img)
            s, rel = base.agreement(ref, a2t_proj(p, cv, dim))
            tried.append((rel, s, cv))
    tried.sort(key=lambda r: r[0])
    rel, s, cv = tried[0]
    print("offset    ASTRA detector shift = {:+d} u, {:+d} v; agree {:.2e} (scale {:.4f}), next best {:.2e}".format(
        cv["off_u"], cv["off_v"], rel, s, tried[1][0]))
    if rel > 0.05 or tried[1][0] < 2 * rel:
        print("  WARNING: the detector offset does not translate clearly: the offset group's 'agree' means nothing.")
    return {"convention": cv, "agree": rel, "scale": s, "next_best_agree": tried[1][0]}


def calibrations_needed(cases):
    beams = {c["beam"] for c in cases}
    return [b for b in ("cone", "parallel", "fan") if b in beams], any(any(c["offdet"]) for c in cases)


def convention_for(c, calib):
    beam = c["beam"]
    if beam == "cone" and any(c["offdet"]):
        return calib["offset"]["convention"]
    return calib[beam]["convention"]


# -------------------------------------------------------------------------------------------------------- the calls

def make_calls(c, mode, devices, conv, niter):
    import tigre
    import tigre.algorithms as algs

    geo = build_geometry(c)
    angles = np.linspace(0, 2 * np.pi, c["angles"], endpoint=False)
    nv, nu = (int(v) for v in geo.nDetector)
    nz, ny, nx = (int(v) for v in geo.nVoxel)
    nang = len(angles)
    vol_geom, proj_geom, dim = astra_setup(c, geo, angles, conv)
    import astra
    astra.set_gpu_index(list(devices))
    A = Astra(vol_geom, proj_geom, dim)
    proj_shape = (nang, nu) if dim == 2 else (nv, nang, nu)
    vol_shape = (ny, nx) if dim == 2 else (nz, ny, nx)
    calls = {}

    if mode == "host":
        ids = base.gpu_ids(devices)
        img = phantom(c["nvox"])
        proj_t = tigre.Ax(img, geo, angles, "interpolated", gpuids=ids)
        empty_proj = lambda: np.empty(proj_shape, dtype=np.float32)  # noqa: E731
        empty_vol = lambda: np.empty(vol_shape, dtype=np.float32)  # noqa: E731
        sync = lambda: None  # noqa: E731
        release = lambda: None  # noqa: E731
        kw = {"gpuids": ids}
    else:
        import torch
        dev = torch.device("cuda", devices[0])
        img = torch.from_numpy(phantom(c["nvox"])).to(dev)
        proj_t = tigre.Ax(img, geo, angles, "interpolated")
        if not (isinstance(proj_t, torch.Tensor) and proj_t.is_cuda):
            raise RuntimeError("this TIGRE does not take GPU tensors (3.3 or later needed)")
        empty_proj = lambda: torch.empty(proj_shape, dtype=torch.float32, device=dev)  # noqa: E731
        empty_vol = lambda: torch.empty(vol_shape, dtype=torch.float32, device=dev)  # noqa: E731
        sync = lambda: torch.cuda.synchronize(dev)  # noqa: E731
        release = torch.cuda.empty_cache  # see tigre_vs_astra.make_calls
        kw = {}
    img_a = t2a_vol(img, dim)
    proj_a = t2a_proj(proj_t, conv, dim)

    def done(r):
        sync()
        return r

    ops = c["ops"]
    if "fp" in ops:
        calls["TIGRE fp Siddon"] = lambda: done(tigre.Ax(img, geo, angles, "Siddon", **kw))
        calls["TIGRE fp interpolated"] = lambda: done(tigre.Ax(img, geo, angles, "interpolated", **kw))

        def astra_fp():
            out = empty_proj()
            A.fp(img_a, out)
            return done(out)
        calls["ASTRA fp"] = astra_fp
    if "bp" in ops:
        if c["beam"] == "parallel":
            calls["TIGRE bp"] = lambda: done(tigre.Atb(proj_t, geo, angles, **kw))
        else:
            calls["TIGRE bp matched"] = lambda: done(tigre.Atb(proj_t, geo, angles, "matched", **kw))
            calls["TIGRE bp FDK"] = lambda: done(tigre.Atb(proj_t, geo, angles, "FDK", **kw))

        def astra_bp():
            out = empty_vol()
            A.bp(out, proj_a)
            return done(out)
        calls["ASTRA bp"] = astra_bp
    if mode == "host":  # ASTRA's algorithms take linked host arrays; TIGRE's take NumPy arrays
        if "fdk" in ops and c["beam"] != "parallel":
            calls["TIGRE fdk"] = lambda: algs.fdk(proj_t, geo, angles, **kw)
            calls["ASTRA fdk"] = lambda: A.algorithm("FBP_CUDA" if dim == 2 else "FDK_CUDA", empty_vol(), proj_a)
        if "sirt" in ops:
            calls["TIGRE sirt"] = lambda: algs.sirt(proj_t, geo, angles, niter, verbose=False, **kw)
            calls["ASTRA sirt"] = lambda: A.algorithm("SIRT_CUDA" if dim == 2 else "SIRT3D_CUDA",
                                                      np.zeros(vol_shape, np.float32), proj_a, niter,
                                                      {"MinConstraint": 0.0})  # TIGRE's noneg=True
        if "cgls" in ops:
            calls["TIGRE cgls"] = lambda: algs.cgls(proj_t, geo, angles, niter, verbose=False, **kw)
            calls["ASTRA cgls"] = lambda: A.algorithm("CGLS_CUDA" if dim == 2 else "CGLS3D_CUDA",
                                                      np.zeros(vol_shape, np.float32), proj_a, niter)

    def checksum():
        return float(base._numpy(proj_a).sum(dtype=np.float64)), float(base._numpy(proj_t).sum(dtype=np.float64))

    return calls, A.close, checksum, release, dim


def bench(c, mode, ngpu, conv, args):
    niter = c["niter"] or 0
    run = {"key": run_key(c, mode, ngpu), "case": c["key"], "group": c["group"], "label": c["label"],
           "beam": c["beam"], "mode": mode, "gpus": ngpu, "nvox": c["nvox"], "ndet": c["ndet"],
           "angles": c["angles"], "niter": niter, "times": {}, "reps": {}, "agree": {}, "errors": {}}
    try:
        calls, cleanup, checksum, release, dim = make_calls(c, mode, list(range(ngpu)), conv, niter)
    except Exception as e:  # e.g. out of memory for the inputs
        run["errors"]["setup"] = "{}: {}".format(type(e).__name__, e)
        print("  setup failed: " + run["errors"]["setup"])
        return run
    try:
        before = checksum()
        ok, first = {}, {}
        for op in OPS:
            aname = "ASTRA " + op
            if aname not in calls:
                continue
            ref, first[aname] = _try_timed(calls[aname], aname, run)
            ok[aname] = ref is not None
            if ref is not None and mode == "device":
                # the check keeps ASTRA's result in host memory: left on the GPU it would take the memory the
                # TIGRE call next to it needs (1024^3: 4.3 GB of 12 GB), and TIGRE would fail where ASTRA ran
                ref = base._numpy(ref)
                release()
            for _, tname, an in PAIRS:
                if an != aname or tname not in calls:
                    continue
                out, first[tname] = _try_timed(calls[tname], tname, run)
                ok[tname] = out is not None
                if out is not None and ref is not None:
                    ref_t = a2t_proj(ref, conv, dim) if op == "fp" else a2t_vol(ref, dim)
                    s, rel = base.agreement(out, ref_t)
                    run["agree"][tname] = {"agree": rel, "scale": s}
                del out
                release()
            del ref
            release()
        names = [k for k, v in ok.items() if v]
        for k in names:  # calls per round: enough for --min-time seconds, at least --reps
            run["reps"][k] = int(min(200, max(args.reps, math.ceil(args.min_time / max(first[k], 1e-6)))))
        run["times"] = {k: [] for k in names}
        for r in range(args.rounds):
            for op in OPS:
                group = [k for k in names if k.split()[1] == op]
                tigre_first = r % 2 == 0
                group.sort(key=lambda k: (k.startswith("ASTRA") == tigre_first, k))
                for k in group:
                    if k not in run["times"]:
                        continue
                    try:
                        run["times"][k].append(base._timed(calls[k], run["reps"][k], release))
                    except Exception as e:
                        run["errors"][k] = "{}: {}".format(type(e).__name__, e)
                        print("  {} failed in round {}: {}".format(k, r + 1, run["errors"][k]))
                        del run["times"][k]
                        release()
            print("    round {}/{}".format(r + 1, args.rounds), flush=True)
        if checksum() != before:
            run["errors"]["inputs"] = "the projections given to bp/fdk/sirt/cgls changed during the run"
    finally:
        cleanup()
    return run


def _try_timed(f, name, run):
    t0 = time.perf_counter()
    try:
        r = f()
    except Exception as e:
        run["errors"][name] = "{}: {}".format(type(e).__name__, e)
        print("  {} failed: {}".format(name, run["errors"][name]))
        return None, None
    return r, time.perf_counter() - t0


# --------------------------------------------------------------------------------------------------------- report

def pair_stats(run, tname, aname):
    return base.pair_stats(run, tname, aname)


def short(tname):
    return tname[6:]


def describe(run):
    nz, ny, nx = run["nvox"]
    nv, nu = run["ndet"]
    img = "{}x{}".format(ny, nx) if run["beam"] == "fan" else "{}x{}x{}".format(nz, ny, nx)
    det = "{}".format(nu) if run["beam"] == "fan" else "{}x{}".format(nv, nu)
    return "{} beam, image {}, detector {}, {} angles{}".format(
        run["beam"], img, det, run["angles"], ", {} iterations".format(run["niter"]) if run["niter"] else "")


def print_report(res):
    m = res["meta"]
    print("\n" + "=" * 110)
    print("TIGRE {} / ASTRA {} | {} | {}".format(m["tigre"], m["astra"], m["platform"], "; ".join(m["gpus"])))
    for name, cal in res.get("calibration", {}).items():
        print("calibration {:<9} agree {:.2e}, scale {:.4f}".format(name, cal["agree"], cal["scale"]))
    groups = list(dict.fromkeys(r["group"] for r in res["runs"]))
    for g in groups:
        runs = [r for r in res["runs"] if r["group"] == g]
        print("\n" + "-" * 110 + "\ngroup: " + g)
        for run in runs:
            print("\n{} | {} arrays, {} GPU{}".format(describe(run), "NumPy" if run["mode"] == "host" else
                                                    "GPU tensor", run["gpus"], "s" if run["gpus"] > 1 else ""))
            print("  {:<18} {:>10} {:>10}  {:<24} {:>16} {:>9} {:>8} {:>6}".format(
                "operation", "TIGRE ms", "ASTRA ms", "TIGRE/ASTRA per round", "mean, 95% CI", "agree", "scale",
                "calls"))
            for _, tname, aname in PAIRS:
                st = pair_stats(run, tname, aname)
                if st is None:
                    continue
                ag = run["agree"].get(tname, {})
                reps = "{}/{}".format(run["reps"].get(tname, "?"), run["reps"].get(aname, "?"))
                print("  {:<18} {:>10.2f} {:>10.2f}  {:<24} {:>7.3f} +- {:<6.3f} {:>9.2e} {:>8.4f} {:>6}".format(
                    short(tname), st["tigre_ms"], st["astra_ms"], " ".join("{:.2f}".format(x) for x in st["ratios"]),
                    st["mean"], st["ci"], ag.get("agree", float("nan")), ag.get("scale", float("nan")), reps))
            for k, v in run["errors"].items():
                print("  ! {}: {}".format(k, v))
        summary(runs)
    scaling(res["runs"])
    print("\nTIGRE/ASTRA < 1: TIGRE faster. agree: relative L2 difference of the outputs after scaling ASTRA's by "
          "'scale' (same problem? not an accuracy measure). calls: calls per round, TIGRE/ASTRA (median kept).")


def summary(runs):
    """One line per TIGRE call, one column per run of the group: the mean ratio."""
    cols = [r for r in runs if r["times"]]
    if len(cols) < 2:
        return
    heads = ["{}{}{}".format(r["label"], ", T" if r["mode"] == "device" else "",
                             ", {}G".format(r["gpus"]) if r["gpus"] > 1 else "") for r in cols]
    width = max(12, max(len(h) for h in heads) + 2)
    print("\n  summary, TIGRE/ASTRA (T: GPU tensors, nG: n GPUs)")
    print("  {:<18}".format("") + "".join("{:>{w}}".format(h, w=width) for h in heads))
    for _, tname, aname in PAIRS:
        cells = [pair_stats(r, tname, aname) for r in cols]
        if not any(cells):
            continue
        print("  {:<18}".format(short(tname)) + "".join(
            "{:>{w}}".format("{:.2f}".format(s["mean"]) if s else "-", w=width) for s in cells))


def scaling(runs):
    """Speed-up from 1 GPU to n GPUs (NumPy arrays), per call, for cases run on several GPU counts."""
    host = [r for r in runs if r["mode"] == "host" and r["times"]]
    by_case = {}
    for r in host:
        by_case.setdefault(r["case"], {})[r["gpus"]] = r
    multi = {k: v for k, v in by_case.items() if 1 in v and len(v) > 1}
    if not multi:
        return
    print("\n" + "-" * 110 + "\nspeed-up from 1 GPU (time on 1 / time on n), NumPy arrays")
    for key, by_g in multi.items():
        counts = sorted(g for g in by_g if g > 1)
        print("\n  " + key)
        print("  {:<24}".format("call") + "".join("{:>10}".format("{} GPUs".format(g)) for g in counts))
        names = [k for k in by_g[1]["times"]]
        for k in sorted(names, key=lambda k: [p[1] for p in PAIRS].index(k) if k.startswith("TIGRE") else
                        [p[2] for p in PAIRS].index(k) + 0.5):
            row = []
            for g in counts:
                a, b = by_g[1]["times"].get(k), by_g[g]["times"].get(k)
                row.append("{:>10.2f}".format(statistics.median(a) / statistics.median(b)) if a and b else
                           "{:>10}".format("-"))
            print("  {:<24}".format(k) + "".join(row))


def write_csv(results, path):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["machine", "group", "label", "beam", "arrays", "gpus", "nz", "ny", "nx", "nv", "nu", "angles",
                    "niter", "tigre_call", "astra_call", "tigre_ms", "astra_ms", "ratio_mean", "ratio_ci95",
                    "agree", "scale"])
        for res in results:
            for r in res["runs"]:
                for _, tname, aname in PAIRS:
                    st = pair_stats(r, tname, aname)
                    if st is None:
                        continue
                    ag = r["agree"].get(tname, {})
                    w.writerow([res["meta"]["label"], r["group"], r["label"], r["beam"],
                                "numpy" if r["mode"] == "host" else "tensor", r["gpus"], *r["nvox"], *r["ndet"],
                                r["angles"], r["niter"], tname, aname, "{:.4f}".format(st["tigre_ms"]),
                                "{:.4f}".format(st["astra_ms"]), "{:.4f}".format(st["mean"]),
                                "{:.4f}".format(st["ci"]), ag.get("agree", ""), ag.get("scale", "")])
    print("wrote " + path)


# ----------------------------------------------------------------------------------------------------------- main

def meta(args, gcounts):
    import astra
    try:
        from importlib.metadata import version
        tigre_version = version("tigre")
    except Exception:
        tigre_version = "?"
    gpus = []
    for i in range(max(gcounts)):
        try:
            gpus.append(str(astra.get_gpu_info(i)).strip())
        except Exception:
            pass
    label = args.label or platform.node().split(".")[0] or "run"
    return {"tigre": tigre_version, "astra": getattr(astra, "__version__", "?"), "platform": platform.platform(),
            "python": sys.version.split()[0], "gpus": gpus or ["?"], "label": re.sub(r"[^A-Za-z0-9]+", "-", label),
            "rounds": args.rounds, "reps": args.reps, "min_time": args.min_time, "groups": args.groups,
            "quick": args.quick, "date": time.strftime("%Y-%m-%d %H:%M")}


def host_bytes(c, mode):
    """A rough upper bound of the host memory one case needs (inputs, outputs of both toolboxes, filtering)."""
    nz, ny, nx = c["nvox"]
    nv, nu = c["ndet"]
    vol, proj = 4 * nz * ny * nx, 4 * nv * nu * c["angles"]
    if mode == "device":
        return vol + 2 * proj
    if set(c["ops"]) <= {"fp", "bp"}:
        # image, projections in TIGRE's and ASTRA's layouts, and the two results compared at a time
        return vol + 2 * proj + 2 * max(vol, proj)
    return 4 * vol + 4 * proj


def total_ram():
    try:
        import psutil
        return psutil.virtual_memory().total
    except Exception:
        pass
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except Exception:
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--groups", nargs="+", default=DEFAULT_GROUPS, choices=GROUPS)
    p.add_argument("--gpus", type=int, nargs="+", default=[1, 2, 4],
                   help="GPU counts for the multi-GPU cases (default 1 2 4; counts above the visible GPUs skipped)")
    p.add_argument("--modes", nargs="+", default=["host", "device"], choices=["host", "device"])
    p.add_argument("--rounds", type=int, default=3)
    p.add_argument("--reps", type=int, default=2, help="least calls per round; the median is kept")
    p.add_argument("--min-time", type=float, default=0.3, help="seconds each measurement lasts at least")
    p.add_argument("--niter", type=int, default=10, help="SIRT/CGLS iterations, cone beam")
    p.add_argument("--fan-niter", type=int, default=50, help="SIRT/CGLS iterations, fan beam")
    p.add_argument("--quick", action="store_true", help="smaller sizes and fewer cases, 2 rounds")
    p.add_argument("--out", default="tigre_vs_astra_suite.json")
    p.add_argument("--resume", action="store_true", help="keep the runs already in --out, run the others")
    p.add_argument("--label", default=None, help="a short name for this machine")
    p.add_argument("--list", action="store_true", help="print the cases and stop")
    p.add_argument("--cases", nargs="+", metavar="KEY", help="only these cases (keys as --list prints them)")
    p.add_argument("--ops", nargs="+", default=None, choices=["fp", "bp", "fdk", "sirt", "cgls"],
                   help="only these operations (default: all; 'fp bp' for the operators alone)")
    p.add_argument("--check", action="store_true", help="only the calibrations")
    p.add_argument("--report", nargs="+", metavar="JSON", help="print saved results instead of running")
    p.add_argument("--csv", metavar="CSV", help="also write every comparison as one CSV row")
    args = p.parse_args()
    if args.quick and "--rounds" not in sys.argv:
        args.rounds = 2

    if args.report:
        results = []
        for f in args.report:
            with open(f) as fh:
                results.append(json.load(fh))
        for res in results:
            print_report(res)
        if args.csv:
            write_csv(results, args.csv)
        return

    cases = build_cases(args.groups, args.quick, args.niter, args.fan_niter)
    if args.cases:
        unknown = set(args.cases) - {c["key"] for c in cases}
        if unknown:
            sys.exit("unknown cases: " + ", ".join(sorted(unknown)))
        cases = [c for c in cases if c["key"] in args.cases]
    if args.ops:
        for c in cases:
            c["ops"] = [o for o in c["ops"] if o in args.ops]
        cases = [c for c in cases if c["ops"]]
    if args.list:
        for c, mode, g in plan(cases, args.gpus, args.modes):
            print("  {:<100} ~{:.1f} GB host".format(run_key(c, mode, g), host_bytes(c, mode) / 2 ** 30))
        return

    if os.path.isdir(os.path.join(os.getcwd(), "Python", "tigre")):
        sys.exit("Run this outside the TIGRE source tree, so that the installed TIGRE is imported.")
    import astra  # noqa: F401
    import tigre  # noqa: F401
    available = base.count_gpus()
    gcounts = sorted({g for g in args.gpus if g <= available} | ({1} if available else set()))
    for g in args.gpus:
        if g > available:
            print("{} GPUs asked, {} visible: skipped".format(g, available))
    if not gcounts:
        sys.exit("no GPU")

    modes = list(args.modes)
    if "device" in modes:
        try:
            import torch
            if not torch.cuda.is_available():
                raise ImportError
        except ImportError:
            print("PyTorch with CUDA not found: GPU-tensor runs skipped")
            modes.remove("device")

    res = {"meta": meta(args, gcounts), "runs": [], "calibration": {}}
    done, old_calibration = set(), {}
    if args.resume and os.path.exists(args.out):
        with open(args.out) as fh:
            old = json.load(fh)
        res["runs"] = [r for r in old["runs"] if "setup" not in r["errors"]]
        old_calibration = old.get("calibration", {})
        done = {r["key"] for r in res["runs"]}
        print("resuming: {} runs kept from {}".format(len(done), args.out))
    print("TIGRE {} / ASTRA {}; {}".format(res["meta"]["tigre"], res["meta"]["astra"], "; ".join(res["meta"]["gpus"])))

    beams, need_offset = calibrations_needed(cases)
    for b in beams:  # always again, also on --resume (a few seconds): a fixed calibration then takes effect
        res["calibration"][b] = calibrate(b)
    if need_offset:
        res["calibration"]["offset"] = calibrate_offset(res["calibration"]["cone"]["convention"])
    if args.resume and done:
        kept_beams = {r["beam"] for r in res["runs"]} | ({"offset"} if any("offset" in k for k in done) else set())
        for b, cal in old_calibration.items():
            if b in kept_beams and b in res["calibration"] and cal["convention"] != res["calibration"][b]["convention"]:
                print("  note: the {} convention differs from the one the kept runs used".format(b))
    if args.check:
        return

    ram = total_ram()
    todo = [(c, m, g) for c, m, g in plan(cases, gcounts, modes) if run_key(c, m, g) not in done]
    t_start = time.time()
    for i, (c, mode, g) in enumerate(todo):
        key = run_key(c, mode, g)
        print("\n[{}/{}] {}".format(i + 1, len(todo), key), flush=True)
        need = host_bytes(c, mode)
        if ram and need > 0.85 * ram:
            print("  skipped: needs about {:.0f} GB of host memory, this machine has {:.0f} GB".format(
                need / 2 ** 30, ram / 2 ** 30))
            continue
        run = bench(c, mode, g, convention_for(c, res["calibration"]), args)
        res["runs"].append(run)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)
        try:
            import torch
            torch.cuda.empty_cache()
        except Exception:
            pass
    order = {run_key(c, m, g): i for i, (c, m, g) in enumerate(plan(cases, gcounts, ["host", "device"]))}
    res["runs"].sort(key=lambda r: order.get(r["key"], len(order)))
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print_report(res)
    if args.csv:
        write_csv([res], args.csv)
    print("\nsaved {} ({:.0f} s)".format(os.path.abspath(args.out), time.time() - t_start))


if __name__ == "__main__":
    main()
