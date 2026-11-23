"""Every comparison of a suite JSON where ASTRA is faster (ratio > 1, and > 1 beyond its 95% interval),
grouped by operation, worst first. Usage: python3 losses.py suite.json [gpus]"""
import json
import statistics as st
import sys

PAIRS = [("TIGRE fp Siddon", "ASTRA fp"), ("TIGRE fp interpolated", "ASTRA fp"), ("TIGRE bp matched", "ASTRA bp"),
         ("TIGRE bp FDK", "ASTRA bp"), ("TIGRE bp", "ASTRA bp"), ("TIGRE fdk", "ASTRA fdk"),
         ("TIGRE sirt", "ASTRA sirt"), ("TIGRE cgls", "ASTRA cgls")]
T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57}


def cause(t, r):
    """First guess at the cause of an ASTRA win, from what is known (to be checked against the profiles)."""
    nz, ny, nx = r["nvox"]
    small = r["beam"] != "fan" and nx <= 128 or r["group"] == "latency" or (r["beam"] == "fan" and nx <= 256)
    if t == "TIGRE fdk":
        return "FDK filtering on the CPU (decision: #423 dropped)"
    if t in ("TIGRE sirt", "TIGRE cgls"):
        return "algorithms in NumPy on the host; CGLS does 2 Ax per iteration (design)"
    if t == "TIGRE fp interpolated":
        return "interpolated kernel: trilinear 3D samples every 0.5 voxel (ignored for now)"
    if small:
        return "fixed per-call cost"
    if t == "TIGRE fp Siddon" and r["beam"] == "fan":
        return "Siddon exact traversal vs ASTRA's Joseph (experiment H)"
    if t == "TIGRE bp" and r["beam"] == "parallel":
        return "parallel bp kernel (experiments F, G)"
    return "unexplained: investigate"


res = json.load(open(sys.argv[1]))
gpus = int(sys.argv[2]) if len(sys.argv) > 2 else 1
rows = []
for r in res["runs"]:
    if r["gpus"] != gpus:
        continue
    for t, a in PAIRS:
        tt, aa = r["times"].get(t), r["times"].get(a)
        if not tt or not aa:
            continue
        q = [x / y for x, y in zip(tt, aa)]
        m = st.mean(q)
        ci = T95.get(len(q) - 1, 2.0) * st.stdev(q) / len(q) ** 0.5 if len(q) > 1 else 0.0
        rows.append((t, r["label"], r["mode"], m, ci, 1000 * st.median(tt), 1000 * st.median(aa), cause(t, r)))
lose = sorted([x for x in rows if x[3] - x[4] > 1], key=lambda x: (x[0], -x[3]))
win = [x for x in rows if x[3] + x[4] < 1]
causes = {}
for x in rows:
    if x[3] - x[4] > 1:
        causes[x[7]] = causes.get(x[7], 0) + 1
tie = [x for x in rows if not (x[3] - x[4] > 1 or x[3] + x[4] < 1)]
print("%d comparisons on %d GPU(s): TIGRE faster %d, ASTRA faster %d, tie %d" % (len(rows), gpus, len(win), len(lose), len(tie)))
cur = None
for t, lab, mode, m, ci, tm, am, why in lose:
    if t != cur:
        print("\n== " + t)
        cur = t
    print("  %5.2f +-%4.2f  %-7s %-30s TIGRE %9.2f ms  ASTRA %9.2f ms  %s" % (m, ci, mode, lab[:30], tm, am, why))
print("\nASTRA wins by cause:")
for k, v in sorted(causes.items(), key=lambda kv: -kv[1]):
    print("  %3d  %s" % (v, k))
