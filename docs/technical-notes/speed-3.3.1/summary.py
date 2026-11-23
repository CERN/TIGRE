"""Summary table for the note: per operation and class of problem, the median TIGRE/ASTRA time ratio (range) of
a suite JSON, and the median speed-up of TIGRE over a baseline suite (3.3.0) on the same cases.
    python3 summary.py suite_rc.json suite_330.json tables/win_summary.tex"""
import collections
import json
import statistics as st
import sys

PAIRS = [("TIGRE fp Siddon", "ASTRA fp", "fp Siddon"), ("TIGRE fp interpolated", "ASTRA fp", "fp interpolated"),
         ("TIGRE bp matched", "ASTRA bp", "bp matched"), ("TIGRE bp FDK", "ASTRA bp", "bp FDK weights"),
         ("TIGRE bp", "ASTRA bp", "bp")]
CLASSES = ["small", "cone", "parallel", "fan"]
NAMES = {"small": "small 3D ($\\le128^3$)", "cone": "cone $256^3$--$1024^3$", "parallel": "parallel",
         "fan": "fan (2D)"}


def load(f):
    d = {}
    for x in json.load(open(f))["runs"]:
        if x["gpus"] != 1:
            continue
        for t, a, _ in PAIRS:
            if t in x["times"] and a in x["times"]:
                d[(x["key"], t)] = (st.median(x["times"][t]), st.median(x["times"][a]), x)
    return d


def cls(x):
    if x["beam"] == "fan":
        return "fan"
    if x["group"] == "latency" or x["nvox"][2] <= 128:
        return "small"
    return x["beam"]


def colour(v):
    s = "%.2f" % v
    return "\\tfast{%s}" % s if v < 0.97 else ("\\tslow{%s}" % s if v > 1.03 else s)


new, old = load(sys.argv[1]), load(sys.argv[2])
g = collections.defaultdict(list)
for k, (t, a, x) in new.items():
    g[(k[1], cls(x), x["mode"])].append((t / a, old[k][0] / t if k in old else None))
rows = []
for t, _, name in PAIRS:
    for c in CLASSES:
        for mode in ("host", "device"):
            v = g.get((t, c, mode))
            if not v:
                continue
            r = [x[0] for x in v]
            o = [x[1] for x in v if x[1]]
            rows.append("%s & %s & %s & %d & %s & %s--%s & %s \\\\" % (
                name, NAMES[c], "NumPy" if mode == "host" else "tensors", len(v), colour(st.median(r)),
                "%.2f" % min(r), "%.2f" % max(r), "%.2f$\\times$" % st.median(o) if o else "--"))
    rows.append("\\addlinespace[2pt]")
with open(sys.argv[3], "w") as fh:
    fh.write("""\\begin{table}[H]
\\centering\\small
\\caption{Summary, one GPU (RTX 4070, Windows): median ratio $t_\\text{TIGRE}/t_\\text{ASTRA}$ over the cases of
each class (below 1: TIGRE faster), its range, and the median speed-up of TIGRE 3.3.1 over 3.3.0 on the same
cases (``--'': not measured with 3.3.0).}
\\label{tab:summary}
\\begin{tabular}{lllrrrr}
\\toprule
operation & problems & arrays & cases & median & range & vs 3.3.0 \\\\
\\midrule
""" + "\n".join(rows[:-1]) + """
\\bottomrule
\\end{tabular}
\\end{table}
""")
print("ok", len(rows))
