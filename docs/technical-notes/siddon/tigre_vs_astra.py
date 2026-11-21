"""TIGRE against the ASTRA Toolbox: speed of the cone-beam operators, on 1 and 2 GPUs.

Runs the installed TIGRE and ASTRA (``pip install astra-toolbox``) in one process, on the same geometry and
the same image, and times:

    fp   forward projection   TIGRE "Siddon" and "interpolated"    against ASTRA FP3D (Joseph-type)
    bp   backprojection       TIGRE "matched" and "FDK" weights     against ASTRA BP3D
    fdk  FDK reconstruction   TIGRE algorithms.fdk                  against ASTRA FDK_CUDA

in two modes:

    host     NumPy arrays in and out (the copies to and from the GPUs are included), on 1 and on 2 GPUs
    device   PyTorch CUDA tensors in and out, no host copies, 1 GPU; fp and bp only (needs PyTorch)

Geometry: TIGRE's default cone beam at n^3 voxels (source-origin 1000 mm, source-detector 1536 mm, a 256 mm
image, an n x n detector of 0.8 mm pixels), as in the Siddon technical note, and 360 angles.

Same problem in both toolboxes. ASTRA's angle and detector conventions differ from TIGRE's. A calibration on a
small problem tries the 8 angle conventions (sign, and an offset of a multiple of 90 degrees) and the 4
detector flips, and keeps the one where ASTRA's projections match TIGRE's. Then, at every size, each TIGRE
output is compared with ASTRA's ("agree": relative L2 difference after a least-squares scale; "scale": that
scale). This only checks that the same problem is timed. It is not an accuracy test: the projection models
differ, and so do the backprojection weights.

Timing. Each round runs every operation of both toolboxes, ``--reps`` times each, keeping the median, and
alternates which toolbox goes first. The report gives the median time over the rounds and, per round, the
ratio TIGRE time / ASTRA time, with its mean and 95% confidence interval (< 1: TIGRE faster).

Usage (run it outside the TIGRE source tree, so that the installed TIGRE is imported):

    python tigre_vs_astra.py --check                          # the calibration only, about a minute
    python tigre_vs_astra.py --out win.json                   # 256^3 and 512^3, on 1 GPU and on 2 if present
    python tigre_vs_astra.py --sizes 256 512 1024 --out linux.json
    python tigre_vs_astra.py --report win.json linux.json --latex astra_results.tex

``--report`` prints saved results again and, with ``--latex``, writes the tables for the technical note.
GPUs are 0 and 1 of the visible ones: use CUDA_VISIBLE_DEVICES to pick others.

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

T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26, 10: 2.23,
       11: 2.20, 12: 2.18, 13: 2.16, 14: 2.14, 15: 2.13, 19: 2.09}

# (operation, TIGRE call, ASTRA call): the pairs that are compared
PAIRS = [
    ("fp", "TIGRE fp Siddon", "ASTRA fp"),
    ("fp", "TIGRE fp interpolated", "ASTRA fp"),
    ("bp", "TIGRE bp matched", "ASTRA bp"),
    ("bp", "TIGRE bp FDK", "ASTRA bp"),
    ("fdk", "TIGRE fdk", "ASTRA fdk"),
]
OPS = ["fp", "bp", "fdk"]


# --------------------------------------------------------------------------------------------- geometry, data

def tigre_geometry(n, ndet=None, ddet=None):
    import tigre
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([n, n, n]))
    if ndet is not None:  # the calibration: a coarse detector that sees the whole image
        geo.nDetector = np.array([ndet, ndet])
        geo.dDetector = np.array([ddet, ddet])
        geo.sDetector = geo.nDetector * geo.dDetector
    return geo


def astra_geometries(geo, angles, conv):
    """ASTRA's volume and projection geometries for TIGRE's geo, with the angle convention conv."""
    import astra
    nz, ny, nx = (int(v) for v in geo.nVoxel)
    sz, sy, sx = (float(v) for v in geo.sVoxel)
    oz, oy, ox = (float(v) for v in geo.offOrigin)
    if any(float(v) != 0 for v in geo.offDetector):
        raise ValueError("detector offsets are not translated to ASTRA")
    vol = astra.create_vol_geom(ny, nx, nz, ox - sx / 2, ox + sx / 2, oy - sy / 2, oy + sy / 2,
                                oz - sz / 2, oz + sz / 2)
    nv, nu = (int(v) for v in geo.nDetector)
    dv, du = (float(v) for v in geo.dDetector)
    a = conv["sign"] * np.asarray(angles, dtype=np.float64) + conv["offset"]
    proj = astra.create_proj_geom("cone", du, dv, nv, nu, a, float(geo.DSO), float(geo.DSD - geo.DSO))
    return vol, proj


def astra_to_tigre(p, conv):
    """ASTRA projections (rows, angles, columns) as TIGRE's (angles, v, u). A view (NumPy or PyTorch)."""
    q = p.transpose(1, 0, 2) if isinstance(p, np.ndarray) else p.permute(1, 0, 2)
    if conv["flip_v"]:
        q = q[:, ::-1, :] if isinstance(q, np.ndarray) else q.flip(1)
    if conv["flip_u"]:
        q = q[:, :, ::-1] if isinstance(q, np.ndarray) else q.flip(2)
    return q


def tigre_to_astra(q, conv):
    """TIGRE projections as ASTRA's, contiguous (ASTRA needs it)."""
    if isinstance(q, np.ndarray):
        if conv["flip_u"]:
            q = q[:, :, ::-1]
        if conv["flip_v"]:
            q = q[:, ::-1, :]
        return np.ascontiguousarray(q.transpose(1, 0, 2))
    if conv["flip_u"]:
        q = q.flip(2)
    if conv["flip_v"]:
        q = q.flip(1)
    return q.permute(1, 0, 2).contiguous()


def phantom(n):
    """A smooth image with no symmetry (Gaussian blobs), so that a wrong orientation shows."""
    c = (np.arange(n, dtype=np.float64) + 0.5) / n - 0.5  # in image widths
    img = np.zeros((n, n, n), dtype=np.float32)
    blobs = [((0.10, -0.20, 0.15), 0.12, 1.0), ((-0.25, 0.10, -0.05), 0.08, 0.7),
             ((0.00, 0.25, 0.30), 0.05, 1.5), ((0.20, 0.00, -0.30), 0.10, 0.5)]
    for (z, y, x), w, amp in blobs:
        gz, gy, gx = (np.exp(-0.5 * ((c - m) / w) ** 2) for m in (z, y, x))
        plane = np.outer(gy, gx).astype(np.float32)
        for k in range(n):  # slice by slice: no n^3 temporaries
            img[k] += np.float32(amp * gz[k]) * plane
    return img


def agreement(t, a):
    """(scale, relative L2 difference of t from scale * a), in double, one slice at a time."""
    t, a = _numpy(t), _numpy(a)
    ta = aa = tt = 0.0
    for i in range(t.shape[0]):
        x = np.asarray(t[i], dtype=np.float64).ravel()
        y = np.asarray(a[i], dtype=np.float64).ravel()
        ta += float(x @ y)
        aa += float(y @ y)
        tt += float(x @ x)
    if aa == 0 or tt == 0:
        return float("nan"), float("nan")
    s = ta / aa
    return s, math.sqrt(max(tt - 2 * s * ta + s * s * aa, 0.0) / tt)


def _numpy(x):
    return x if isinstance(x, np.ndarray) else x.detach().cpu().numpy()


def gpu_ids(devices):
    from tigre.utilities.gpu import GpuIds
    ids = GpuIds()
    ids.devices = list(devices)
    return ids


def count_gpus():
    try:
        import torch
        if torch.cuda.is_available():
            return torch.cuda.device_count()
    except ImportError:
        pass
    try:
        from tigre.utilities.gpu import GpuIds
        return len(GpuIds().devices)
    except Exception:
        return 0


# ------------------------------------------------------------------------------------------------ calibration

def calibrate():
    """The ASTRA angle convention and detector flips under which ASTRA's projections match TIGRE's."""
    import astra
    import astra.experimental as ae
    import tigre
    n, nang = 64, 24
    geo = tigre_geometry(n, ndet=64, ddet=6.4)
    angles = np.linspace(0, 2 * np.pi, nang, endpoint=False)
    img = phantom(n)
    ref = tigre.Ax(img, geo, angles, "interpolated", gpuids=gpu_ids([0]))
    astra.set_gpu_index([0])
    tried = []
    for sign in (1, -1):
        for k in range(4):
            base = {"sign": sign, "offset": k * math.pi / 2}
            vol, proj = astra_geometries(geo, angles, dict(base, flip_u=False, flip_v=False))
            pid = astra.create_projector("cuda3d", proj, vol)
            try:
                p = np.empty((n, nang, n), dtype=np.float32)
                ae.direct_FP3D(pid, img, p)
            finally:
                astra.projector3d.delete(pid)
            for flip_v in (False, True):
                for flip_u in (False, True):
                    conv = dict(base, flip_u=flip_u, flip_v=flip_v)
                    s, rel = agreement(ref, astra_to_tigre(p, conv))
                    tried.append((rel, s, conv))
    tried.sort(key=lambda r: r[0])
    rel, s, conv = tried[0]
    out = {"convention": conv, "agree": rel, "scale": s, "next_best_agree": tried[1][0]}
    print("calibration: ASTRA angle = {:+d} * TIGRE angle + {:.0f} deg, flip u: {}, flip v: {}".format(
        conv["sign"], math.degrees(conv["offset"]), conv["flip_u"], conv["flip_v"]))
    print("  agree {:.2e} (scale {:.4f}); next best convention {:.2e}".format(rel, s, tried[1][0]))
    if rel > 0.05 or tried[1][0] < 2 * rel:
        print("  WARNING: no convention clearly matches. The timings still stand, but the 'agree' columns "
              "do not mean the same problem was solved.")
    return out


# --------------------------------------------------------------------------------------------------- the calls

def make_calls(mode, n, nang, devices, conv, ops):
    """{name: function} for one size and set of GPUs, plus the data they hold."""
    import astra
    import astra.experimental as ae
    import tigre
    import tigre.algorithms as algs

    geo = tigre_geometry(n)
    angles = np.linspace(0, 2 * np.pi, nang, endpoint=False)
    nv, nu = (int(v) for v in geo.nDetector)
    vol_geom, proj_geom = astra_geometries(geo, angles, conv)
    astra.set_gpu_index(list(devices))
    pid = astra.create_projector("cuda3d", proj_geom, vol_geom)
    calls = {}

    if mode == "host":
        ids = gpu_ids(devices)
        img = phantom(n)
        proj_t = tigre.Ax(img, geo, angles, "interpolated", gpuids=ids)  # input of bp and fdk
        proj_a = tigre_to_astra(proj_t, conv)

        def empty_proj():
            return np.empty((nv, nang, nu), dtype=np.float32)

        def empty_vol():
            return np.empty((n, n, n), dtype=np.float32)

        sync = lambda: None  # noqa: E731
        release = lambda: None  # noqa: E731
        kw = {"gpuids": ids}
    else:
        import torch
        dev = torch.device("cuda", devices[0])
        img = torch.from_numpy(phantom(n)).to(dev)
        proj_t = tigre.Ax(img, geo, angles, "interpolated")
        if not (isinstance(proj_t, torch.Tensor) and proj_t.is_cuda):
            raise RuntimeError("this TIGRE does not take GPU tensors (3.3 or later needed)")
        proj_a = tigre_to_astra(proj_t, conv)

        def empty_proj():
            return torch.empty((nv, nang, nu), dtype=torch.float32, device=dev)

        def empty_vol():
            return torch.empty((n, n, n), dtype=torch.float32, device=dev)

        sync = lambda: torch.cuda.synchronize(dev)  # noqa: E731
        # Outputs freed by Python go back to PyTorch's cache, where TIGRE and ASTRA (cudaMalloc) cannot use them:
        # return them to the GPU between calls (not timed).
        release = torch.cuda.empty_cache
        kw = {}

    def done(r):
        sync()
        return r

    if "fp" in ops:
        calls["TIGRE fp Siddon"] = lambda: done(tigre.Ax(img, geo, angles, "Siddon", **kw))
        calls["TIGRE fp interpolated"] = lambda: done(tigre.Ax(img, geo, angles, "interpolated", **kw))

        def astra_fp():
            out = empty_proj()
            ae.direct_FP3D(pid, img, out)
            return done(out)
        calls["ASTRA fp"] = astra_fp
    if "bp" in ops:
        calls["TIGRE bp matched"] = lambda: done(tigre.Atb(proj_t, geo, angles, "matched", **kw))
        calls["TIGRE bp FDK"] = lambda: done(tigre.Atb(proj_t, geo, angles, "FDK", **kw))

        def astra_bp():
            out = empty_vol()
            ae.direct_BP3D(pid, out, proj_a)
            return done(out)
        calls["ASTRA bp"] = astra_bp
    if "fdk" in ops and mode == "host":
        calls["TIGRE fdk"] = lambda: algs.fdk(proj_t, geo, angles, **kw)

        def astra_fdk():
            out = empty_vol()
            vid = astra.data3d.link("-vol", vol_geom, out)
            sid = astra.data3d.link("-sino", proj_geom, proj_a)
            cfg = astra.astra_dict("FDK_CUDA")
            cfg["ReconstructionDataId"] = vid
            cfg["ProjectionDataId"] = sid
            aid = astra.algorithm.create(cfg)
            try:
                astra.algorithm.run(aid)
            finally:
                astra.algorithm.delete(aid)
                astra.data3d.delete([vid, sid])
            return out
        calls["ASTRA fdk"] = astra_fdk

    def cleanup():
        astra.projector3d.delete(pid)

    def input_checksum():
        return float(_numpy(proj_a).sum(dtype=np.float64)), float(_numpy(proj_t).sum(dtype=np.float64))

    return calls, cleanup, input_checksum, conv, release


def bench(mode, n, nang, devices, conv, ops, rounds, reps):
    run = {"mode": mode, "n": n, "angles": nang, "gpus": len(devices), "times": {}, "agree": {}, "errors": {}}
    try:
        calls, cleanup, checksum, conv, release = make_calls(mode, n, nang, devices, conv, ops)
    except Exception as e:  # e.g. out of memory for the inputs
        run["errors"]["setup"] = "{}: {}".format(type(e).__name__, e)
        print("  setup failed: " + run["errors"]["setup"])
        return run
    try:
        before = checksum()
        # Warm-up, one call of each, and the comparison of the outputs: ASTRA's first, then each TIGRE one.
        ok = {}
        for op in ops:
            astra_name = "ASTRA " + op
            if astra_name not in calls:
                continue
            ref = _try(calls[astra_name], astra_name, run)
            ok[astra_name] = ref is not None
            for _, tname, aname in PAIRS:
                if aname != astra_name or tname not in calls:
                    continue
                out = _try(calls[tname], tname, run)
                ok[tname] = out is not None
                if out is not None and ref is not None:
                    if op == "fp":
                        ref_t = astra_to_tigre(ref, conv)
                    else:
                        ref_t = ref
                    s, rel = agreement(out, ref_t)
                    run["agree"][tname] = {"agree": rel, "scale": s}
                del out
                release()
            del ref
            release()
        names = [k for k, v in ok.items() if v]
        run["times"] = {k: [] for k in names}
        for r in range(rounds):
            for op in ops:
                group = [k for k in names if k.split()[1] == op]
                tigre_first = r % 2 == 0
                group.sort(key=lambda k: (k.startswith("ASTRA") == tigre_first, k))
                for k in group:
                    if k not in run["times"]:
                        continue  # failed in an earlier round
                    try:
                        run["times"][k].append(_timed(calls[k], reps, release))
                    except Exception as e:  # e.g. out of GPU memory: drop this call, keep the others
                        run["errors"][k] = "{}: {}".format(type(e).__name__, e)
                        print("  {} failed in round {}: {}".format(k, r + 1, run["errors"][k]))
                        del run["times"][k]
                        release()
            print("    round {}/{}".format(r + 1, rounds), flush=True)
        after = checksum()
        if after != before:
            run["errors"]["inputs"] = "the projections given to bp/fdk changed during the run: {} -> {}".format(
                before, after)
    finally:
        cleanup()
    return run


def _try(f, name, run):
    try:
        return f()
    except Exception as e:
        run["errors"][name] = "{}: {}".format(type(e).__name__, e)
        print("  {} failed: {}".format(name, run["errors"][name]))
        return None


def _timed(f, reps, release):
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        r = f()
        ts.append(time.perf_counter() - t0)
        del r
        release()
    return statistics.median(ts)


# ------------------------------------------------------------------------------------------------------ report

def ci(values):
    m = statistics.mean(values)
    if len(values) < 2:
        return m, float("nan")
    return m, T95.get(len(values) - 1, 2.0) * statistics.stdev(values) / math.sqrt(len(values))


def pair_stats(run, tname, aname):
    t, a = run["times"].get(tname), run["times"].get(aname)
    if not t or not a:
        return None
    ratios = [x / y for x, y in zip(t, a)]
    mean, half = ci(ratios)
    return {"tigre_ms": 1000 * statistics.median(t), "astra_ms": 1000 * statistics.median(a),
            "ratios": ratios, "mean": mean, "ci": half}


def print_report(res):
    m = res["meta"]
    print("\n" + "=" * 100)
    print("TIGRE {} / ASTRA {} | {} | {}".format(m["tigre"], m["astra"], m["platform"], "; ".join(m["gpus"])))
    if res.get("calibration"):
        c = res["calibration"]
        print("calibration agree {:.2e}, scale {:.4f}".format(c["agree"], c["scale"]))
    for run in res["runs"]:
        head = "{} arrays, {} GPU{}, {}^3, {} angles".format(
            "NumPy" if run["mode"] == "host" else "GPU tensor", run["gpus"], "s" if run["gpus"] > 1 else "",
            run["n"], run["angles"])
        print("\n" + head)
        print("  {:<22} {:>10} {:>10}  {:<30} {:>16} {:>9} {:>8}".format(
            "operation", "TIGRE ms", "ASTRA ms", "TIGRE/ASTRA each round", "mean, 95% CI", "agree", "scale"))
        for _, tname, aname in PAIRS:
            st = pair_stats(run, tname, aname)
            if st is None:
                continue
            ag = run["agree"].get(tname, {})
            print("  {:<22} {:>10.1f} {:>10.1f}  {:<30} {:>7.3f} +- {:<6.3f} {:>9.2e} {:>8.4f}".format(
                tname[6:], st["tigre_ms"], st["astra_ms"], " ".join("{:.2f}".format(x) for x in st["ratios"]),
                st["mean"], st["ci"], ag.get("agree", float("nan")), ag.get("scale", float("nan"))))
        for k, v in run["errors"].items():
            print("  ! {}: {}".format(k, v))
    # Multi-GPU scaling: time on 1 GPU / time on 2 GPUs, per call
    host = [r for r in res["runs"] if r["mode"] == "host"]
    sizes = sorted({r["n"] for r in host})
    one = {r["n"]: r for r in host if r["gpus"] == 1}
    two = {r["n"]: r for r in host if r["gpus"] == 2}
    if any(n in one and n in two for n in sizes):
        print("\nspeed-up from 1 to 2 GPUs (time on 1 / time on 2), NumPy arrays")
        names = list(dict.fromkeys(x for _, tn, an in PAIRS for x in (tn, an)))
        print("  {:<22}".format("call") + "".join("{:>10}".format("{}^3".format(n)) for n in sizes))
        for k in names:
            row = []
            for n in sizes:
                a, b = one.get(n, {}).get("times", {}).get(k), two.get(n, {}).get("times", {}).get(k)
                row.append("{:>10.2f}".format(statistics.median(a) / statistics.median(b)) if a and b else
                           "{:>10}".format("-"))
            print("  {:<22}".format(k) + "".join(row))
    print("\nTIGRE/ASTRA < 1: TIGRE faster. agree: relative L2 difference between the two outputs after "
          "scaling ASTRA's by 'scale' (a check that the same problem was solved, not an accuracy measure).")


def latex_tables(results, ops):
    """One table per results file: the pairs of the operations in ops, NumPy rows then GPU-tensor rows."""
    out = []
    labels = {"TIGRE fp Siddon": "fp, Siddon", "TIGRE fp interpolated": "fp, interpolated",
              "TIGRE bp matched": "bp, matched", "TIGRE bp FDK": "bp, FDK weights", "TIGRE fdk": "FDK"}
    for res in results:
        m = res["meta"]
        gcounts = sorted({r["gpus"] for r in res["runs"] if r["mode"] == "host"})
        cols = "ll" + "rrc" * len(gcounts)
        lines = [r"\begin{table}[h]", r"\centering",
                 r"\caption{\label{tab:astra-%s}TIGRE against ASTRA, %s (%s). Time of the whole call, %d angles; "
                 r"median over %d rounds. Ratio: TIGRE time over ASTRA time, mean and 95\%% confidence interval "
                 r"over rounds (below 1: TIGRE faster). ``tensors'': PyTorch GPU tensors in and out, no host "
                 r"copies; otherwise NumPy arrays.}" % (
                     m["label"], _tex(gpu_name(m["gpus"][0])), _tex(m["platform"].split("-")[0]),
                     res["runs"][0]["angles"] if res["runs"] else 0, m["rounds"]),
                 r"\small\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{%s}" % cols, r"\toprule"]
        lines.append("& & " + " & ".join(r"\multicolumn{3}{c}{%d GPU%s}" % (g, "s" if g > 1 else "")
                                         for g in gcounts) + r"\\")
        lines.append("".join(r"\cmidrule(lr){%d-%d}" % (3 + 3 * i, 5 + 3 * i) for i in range(len(gcounts))))
        lines.append("image & TIGRE operation & " + " & ".join(["TIGRE ms & ASTRA ms & ratio"] * len(gcounts))
                     + r"\\")
        lines.append(r"\midrule")
        for mode in ("host", "device"):
            for n in sorted({r["n"] for r in res["runs"] if r["mode"] == mode}):
                first = True
                for op, tname, aname in PAIRS:
                    if op not in ops:
                        continue
                    cells = []
                    any_value = False
                    for g in gcounts:
                        run = next((r for r in res["runs"] if r["mode"] == mode and r["n"] == n
                                    and r["gpus"] == g), None)
                        st = pair_stats(run, tname, aname) if run else None
                        if st is None:
                            cells.append("--- & --- & ---")
                        else:
                            any_value = True
                            cells.append(r"%.0f & %.0f & $%.2f\pm%.2f$" % (st["tigre_ms"], st["astra_ms"],
                                                                         st["mean"], st["ci"]))
                    if not any_value:
                        continue
                    image = (r"$%d^3$%s" % (n, ", tensors" if mode == "device" else "")) if first else ""
                    first = False
                    lines.append("%s & %s & %s\\\\" % (image, labels[tname], " & ".join(cells)))
        lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
        out.append("\n".join(lines))
    return "\n".join(out)


def gpu_name(info):
    """'GPU #0: NVIDIA GeForce RTX 4070, with ...' -> 'NVIDIA GeForce RTX 4070'."""
    return re.sub(r"^GPU\s*#?\d+:\s*", "", info.split(",")[0]).strip()


def _tex(s):
    return s.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&").replace("#", r"\#")


# -------------------------------------------------------------------------------------------------------- main

def meta(args):
    import astra
    try:
        from importlib.metadata import version
        tigre_version = version("tigre")
    except Exception:
        tigre_version = "?"
    gpus = []
    for i in range(max(args.gpus)):
        try:
            gpus.append(str(astra.get_gpu_info(i)).strip())
        except Exception:
            pass
    label = args.label or platform.node().split(".")[0] or "run"
    return {"tigre": tigre_version, "astra": getattr(astra, "__version__", "?"), "platform": platform.platform(),
            "python": sys.version.split()[0], "gpus": gpus or ["?"], "label": re.sub(r"[^A-Za-z0-9]+", "-", label),
            "rounds": args.rounds, "reps": args.reps, "date": time.strftime("%Y-%m-%d %H:%M")}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--sizes", type=int, nargs="+", default=[256, 512])
    p.add_argument("--angles", type=int, default=360)
    p.add_argument("--gpus", type=int, nargs="+", default=[1, 2], help="GPU counts to run (default 1 2)")
    p.add_argument("--ops", nargs="+", default=OPS, choices=OPS)
    p.add_argument("--modes", nargs="+", default=["host", "device"], choices=["host", "device"])
    p.add_argument("--rounds", type=int, default=4)
    p.add_argument("--reps", type=int, default=2, help="calls per round; the median is kept")
    p.add_argument("--out", default="tigre_vs_astra.json")
    p.add_argument("--label", default=None, help="a short name for this machine, used in table labels")
    p.add_argument("--check", action="store_true", help="only the calibration")
    p.add_argument("--report", nargs="+", metavar="JSON", help="print saved results instead of running")
    p.add_argument("--latex", metavar="TEX", help="with --report: write the LaTeX tables here")
    p.add_argument("--latex-ops", nargs="+", default=["fp"], choices=OPS)
    args = p.parse_args()

    if args.report:
        results = []
        for f in args.report:
            with open(f) as fh:
                results.append(json.load(fh))
        for res in results:
            print_report(res)
        if args.latex:
            with open(args.latex, "w") as fh:
                fh.write(latex_tables(results, args.latex_ops))
            print("\nwrote " + args.latex)
        return

    if os.path.isdir(os.path.join(os.getcwd(), "Python", "tigre")):
        sys.exit("Run this outside the TIGRE source tree, so that the installed TIGRE is imported.")
    import astra  # noqa: F401  (fail early if missing)
    import tigre  # noqa: F401
    available = count_gpus()
    gcounts = [g for g in args.gpus if g <= available]
    for g in args.gpus:
        if g > available:
            print("{} GPUs asked, {} visible: skipped".format(g, available))
    if not gcounts:
        sys.exit("no GPU")
    args.gpus = gcounts
    res = {"meta": meta(args), "runs": []}
    print("TIGRE {} / ASTRA {}; {}".format(res["meta"]["tigre"], res["meta"]["astra"], "; ".join(res["meta"]["gpus"])))

    res["calibration"] = calibrate()
    conv = res["calibration"]["convention"]
    if args.check:
        return

    modes = list(args.modes)
    if "device" in modes:
        try:
            import torch
            if not torch.cuda.is_available():
                raise ImportError
        except ImportError:
            print("PyTorch with CUDA not found: GPU-tensor mode skipped")
            modes.remove("device")
    plan = [("host", n, g) for n in args.sizes for g in gcounts]
    if "device" in modes:
        plan += [("device", n, 1) for n in args.sizes]
    if "host" not in modes:
        plan = [x for x in plan if x[0] != "host"]
    for mode, n, g in plan:
        print("\n{} arrays, {} GPU(s), {}^3, {} angles".format("NumPy" if mode == "host" else "GPU tensor", g, n,
                                                           args.angles), flush=True)
        res["runs"].append(bench(mode, n, args.angles, list(range(g)), conv, args.ops, args.rounds, args.reps))
        with open(args.out, "w") as fh:  # after every run, so that a crash keeps what was measured
            json.dump(res, fh, indent=1)
        try:
            import torch
            torch.cuda.empty_cache()
        except Exception:
            pass
    print_report(res)
    print("\nsaved " + os.path.abspath(args.out))


if __name__ == "__main__":
    main()
