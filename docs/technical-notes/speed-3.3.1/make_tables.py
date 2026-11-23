"""LaTeX tables for the TIGRE vs ASTRA suite report, from suite.json (one file per table, in tables/)."""
import json
import os
import statistics as st
import sys

PAIRS = [
    ("TIGRE fp Siddon", "ASTRA fp", "fp, Siddon"),
    ("TIGRE fp interpolated", "ASTRA fp", "fp, interpolated"),
    ("TIGRE bp matched", "ASTRA bp", "bp, matched"),
    ("TIGRE bp FDK", "ASTRA bp", "bp, FDK weights"),
    ("TIGRE bp", "ASTRA bp", "bp"),
]
T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57}

res = json.load(open(sys.argv[1]))
out = sys.argv[2]
PREFIX = sys.argv[3] if len(sys.argv) > 3 else ""
MULTI = any(r["gpus"] > 1 for r in json.load(open(sys.argv[1]))["runs"])
os.makedirs(out, exist_ok=True)
runs = res["runs"]


def find(case_prefix, mode, g):
    for r in runs:
        if r["case"] == case_prefix and r["mode"] == mode and r["gpus"] == g:
            return r
    return None


def stats(r, t, a):
    if r is None or not r["times"].get(t) or not r["times"].get(a):
        return None
    tt, aa = r["times"][t], r["times"][a]
    ratios = [x / y for x, y in zip(tt, aa)]
    m = st.mean(ratios)
    ci = T95.get(len(ratios) - 1, 2.0) * st.stdev(ratios) / len(ratios) ** 0.5 if len(ratios) > 1 else float("nan")
    return {"t": 1000 * st.median(tt), "a": 1000 * st.median(aa), "m": m, "ci": ci}


def ms(x):
    if x < 10:
        return "%.2f" % x
    if x < 100:
        return "%.1f" % x
    return "%.0f" % x


def ratio(s):
    if s is None:
        return "---"
    mark = r"$^\dagger$" if s["m"] - s["ci"] <= 1 <= s["m"] + s["ci"] else ""
    cmd = r"\tfast" if s["m"] < 1 else r"\tslow"
    ci = ("%.2f" % s["ci"]) if s["ci"] < 10 else ("%.0f" % s["ci"])
    return r"%s{$%.2f\pm%s$}%s" % (cmd, s["m"], ci, mark)


def cases_of(group):
    seen = []
    for r in runs:
        if r["group"] == group and r["case"] not in seen:
            seen.append(r["case"])
    return seen


def label_of(case):
    return next(r["label"] for r in runs if r["case"] == case)


def tex(s):
    return s.replace("^2", r"$^2$").replace("^3", r"$^3$").replace("_", r"\_")


def one_gpu_table(group, name, caption, with_tensors=True, sep=4):
    cols = "llrrr" + ("rrr" if with_tensors else "")
    L = [r"\begin{table}[H]", r"\centering", r"\caption{%s}" % caption, r"\label{tab:%s}" % name,
         r"\small\setlength{\tabcolsep}{%gpt}" % sep, r"\begin{adjustbox}{max width=\linewidth}", r"\begin{tabular}{%s}" % cols, r"\toprule"]
    if with_tensors:
        L.append(r"& & \multicolumn{3}{c}{NumPy arrays, 1 GPU} & \multicolumn{3}{c}{PyTorch GPU tensors, 1 GPU}\\")
        L.append(r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}")
        L.append(r"problem & operation & TIGRE ms & ASTRA ms & ratio & TIGRE ms & ASTRA ms & ratio\\")
    else:
        L.append(r"& & \multicolumn{3}{c}{NumPy arrays, 1 GPU}\\")
        L.append(r"\cmidrule(lr){3-5}")
        L.append(r"problem & operation & TIGRE ms & ASTRA ms & ratio\\")
    L.append(r"\midrule")
    first_case = True
    for case in cases_of(group):
        h, d = find(case, "host", 1), find(case, "device", 1)
        rows = []
        for t, a, lab in PAIRS:
            sh, sd = stats(h, t, a), stats(d, t, a)
            if sh is None and sd is None:
                continue
            cells = [ms(sh["t"]), ms(sh["a"]), ratio(sh)] if sh else ["---"] * 3
            if with_tensors:
                cells += [ms(sd["t"]), ms(sd["a"]), ratio(sd)] if sd else ["---"] * 3
            rows.append((lab, cells))
        if not rows:
            continue
        if not first_case:
            L.append(r"\addlinespace")
        first_case = False
        for i, (lab, cells) in enumerate(rows):
            L.append("%s & %s & %s\\\\" % (tex(label_of(case)) if i == 0 else "", lab, " & ".join(cells)))
    L += [r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}", r"\end{table}", ""]
    open(os.path.join(out, PREFIX + name + ".tex"), "w").write("\n".join(L))


def speedup(case, k, g):
    a, b = find(case, "host", 1), find(case, "host", g)
    if a is None or b is None or not a["times"].get(k) or not b["times"].get(k):
        return None
    return st.median(a["times"][k]) / st.median(b["times"][k])


def multi_gpu_table(group, name, caption, cases=None):
    L = [r"\begin{table}[H]", r"\centering", r"\caption{%s}" % caption, r"\label{tab:%s}" % name,
         r"\small\setlength{\tabcolsep}{3.5pt}", r"\begin{adjustbox}{max width=\linewidth}", r"\begin{tabular}{llrrrrrrr}", r"\toprule",
         r"& & \multicolumn{3}{c}{TIGRE/ASTRA ratio} & \multicolumn{2}{c}{TIGRE speed-up} & "
         r"\multicolumn{2}{c}{ASTRA speed-up}\\",
         r"\cmidrule(lr){3-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}",
         r"problem & operation & 1 GPU & 2 GPUs & 4 GPUs & 2 GPUs & 4 GPUs & 2 GPUs & 4 GPUs\\", r"\midrule"]
    first_case = True
    for case in (cases or cases_of(group)):
        rows = []
        for t, a, lab in PAIRS:
            s = [stats(find(case, "host", g), t, a) for g in (1, 2, 4)]
            if s[0] is None:
                continue
            sp = [speedup(case, t, 2), speedup(case, t, 4), speedup(case, a, 2), speedup(case, a, 4)]
            rows.append((lab, [ratio(x) for x in s] + ["%.2f" % x if x else "---" for x in sp]))
        if not rows:
            continue
        if not first_case:
            L.append(r"\addlinespace")
        first_case = False
        for i, (lab, cells) in enumerate(rows):
            L.append("%s & %s & %s\\\\" % (tex(label_of(case)) if i == 0 else "", lab, " & ".join(cells)))
    L += [r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}", r"\end{table}", ""]
    open(os.path.join(out, PREFIX + name + ".tex"), "w").write("\n".join(L))


def fixed_cost_table(name, caption, items):
    """t(g) = fixed + parallel/g fitted to the 1- and 4-GPU times; the 2-GPU time checks the fit."""
    L = [r"\begin{table}[H]", r"\centering", r"\caption{%s}" % caption, r"\label{tab:%s}" % name,
         r"\small\setlength{\tabcolsep}{4pt}", r"\begin{adjustbox}{max width=\linewidth}", r"\begin{tabular}{llrrrrrr}", r"\toprule",
         r"problem & call & $t_1$ ms & $t_4$ ms & fixed ms & fixed \% & $t_2$ ms & $t_2$ fit ms\\", r"\midrule"]
    first = True
    for case, calls in items:
        if not first:
            L.append(r"\addlinespace")
        first = False
        for i, k in enumerate(calls):
            t = {g: st.median(find(case, "host", g)["times"][k]) * 1000 for g in (1, 2, 4)}
            p = (t[1] - t[4]) / 0.75
            s = t[1] - p
            L.append("%s & %s & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f\\\\" % (
                (case.split(":")[1].split()[0] + " " + tex(label_of(case))) if i == 0 else "", k, t[1], t[4], s, 100 * s / t[1], t[2], s + p / 2))
    L += [r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}", r"\end{table}", ""]
    open(os.path.join(out, PREFIX + name + ".tex"), "w").write("\n".join(L))


def calib_table():
    L = [r"\begin{table}[H]", r"\centering",
         r"\caption{Calibration: the ASTRA convention that matches TIGRE, found on a small problem for each beam "
         r"type, and how well the projections agree under it (relative $L_2$ difference after a least-squares "
         r"scale). ``Next best'': the best convention that gives measurably different projections.}",
         r"\label{tab:calib}", r"\small", r"\begin{adjustbox}{max width=\linewidth}", r"\begin{tabular}{llrrr}", r"\toprule",
         r"beam & ASTRA convention & agree & scale & next best\\", r"\midrule"]
    for b, c in res["calibration"].items():
        v = c["convention"]
        if b == "offset":
            conv = r"detector shift $%+d\,u, %+d\,v$" % (v["off_u"], v["off_v"])
        else:
            conv = r"$\theta_A = %+d\,\theta_T %+d^\circ$%s" % (
                v["sign"], round(v["offset"] * 180 / 3.141592653589793), ", $u$ flipped" if v["flip_u"] else "")
        nb = c["next_best_agree"]
        nb_s = r"%.1e (tie, see text)" % nb if nb < 1.01 * c["agree"] + 1e-9 else "%.1e" % nb
        L.append(r"%s & %s & %.1e & %.4f & %s\\" % (b, conv, c["agree"], c["scale"], nb_s))
    L += [r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}", r"\end{table}", ""]
    open(os.path.join(out, PREFIX + "calib.tex"), "w").write("\n".join(L))


one_gpu_table("latency", "latency", "Small problems (cone beam, $n^3$ image, $n^2$ detector, 100 angles), 1 GPU.")
one_gpu_table("fan", "fan", "Fan beam (2D): TIGRE's one-slice cone geometry against ASTRA's \\texttt{fanflat} "
              "(FP\\_CUDA, BP\\_CUDA, FBP\\_CUDA), $n^2$ image, $n$ detector pixels, 1 GPU.", sep=3)
one_gpu_table("size", "size1", "Cone beam, $n^3$ image, $n^2$ detector covering the image (409.6\\,mm), "
              "360 angles, 1 GPU.")
MULTI and multi_gpu_table("size", "size_multi", "Cone beam size sweep on 1, 2 and 4 GPUs (NumPy arrays). Speed-up: time on "
                "1 GPU over time on $g$ GPUs, for each toolbox on its own.")
one_gpu_table("legacy", "legacy", "The 3.3.0 note's detector (0.8\\,mm pixels whatever $n$), 1 GPU. Compare with "
              "Table~\\ref{tab:size1} at the same $n$.")
one_gpu_table("angles", "angles", "Cone beam $512^3$, varying the number of angles, 1 GPU.")
one_gpu_table("shape", "shape", "Non-cubic problems, 1 GPU: clinical CBCT ($256\\times512\\times512$ image, "
              "$512\\times1024$ detector), fine detector ($256^3$ image, $1024^2$ detector), coarse detector "
              "($512^3$ image, $256^2$ detector), 360 angles.")
one_gpu_table("offset", "offset", "Cone beam $512^3$ with the detector shifted sideways by 102.4\\,mm "
              "(a quarter of its width; ASTRA \\texttt{cone\\_vec}), 1 GPU.")
one_gpu_table("parallel", "parallel1", "Parallel beam (TIGRE parallel mode, ASTRA \\texttt{parallel3d}), $n^3$, "
              "360 angles, 1 GPU.")
MULTI and multi_gpu_table("parallel", "parallel_multi", "Parallel beam on 1, 2 and 4 GPUs (NumPy arrays).")
# operators only: the reconstructions (FDK, SIRT, CGLS) are not part of the comparison
MULTI and fixed_cost_table("fixed", "Each call's time split into a part that does not shrink with more GPUs (``fixed'') and a "
                 "part that does, fitted as $t_g = \\text{fixed} + \\text{parallel}/g$ to the 1- and 4-GPU times "
                 "(NumPy arrays). The measured 2-GPU time $t_2$ checks the fit.",
                 [(c, k) for c, k in [
                     ("size: cone 512x512x512, det 512x512, 360 angles",
                      ["TIGRE fp Siddon", "TIGRE fp interpolated", "TIGRE bp matched", "ASTRA fp", "ASTRA bp"]),
                     ("size: cone 1024x1024x1024, det 1024x1024, 360 angles",
                      ["TIGRE fp Siddon", "TIGRE fp interpolated", "TIGRE bp matched", "TIGRE bp FDK",
                       "ASTRA fp", "ASTRA bp"]),
                     ("parallel: parallel 1024x1024x1024, det 1024x1024, 360 angles",
                      ["TIGRE fp Siddon", "TIGRE fp interpolated", "TIGRE bp", "ASTRA fp", "ASTRA bp"])]])
calib_table()
print("ok")
