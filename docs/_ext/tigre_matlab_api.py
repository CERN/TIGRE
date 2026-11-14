"""MATLAB API reference pages, from the help text of TIGRE's .m files.

For each group in ``GROUPS`` this writes ``docs/_generated/matlab/<group>.md``
(a MyST fragment), which ``docs/api/matlab/<group>.md`` includes. Each function
gets a heading, its signature, its help text (what ``help NAME`` prints in
MATLAB, without the license block) and a link to the source.
"""
from __future__ import annotations

import pathlib
import re

DOCS = pathlib.Path(__file__).resolve().parent.parent
MATLAB = DOCS.parent / "MATLAB"
OUT = DOCS / "_generated" / "matlab"
SOURCE_URL = "https://github.com/CERN/TIGRE/blob/master/MATLAB/"

# group -> list of paths relative to MATLAB/ (a folder means every .m file in it)
GROUPS = {
    "operators": ["Utilities/Ax.m", "Utilities/Atb.m"],
    "algorithms": ["Algorithms"],
    "geometry": [
        "Utilities/defaultGeometry.m",
        "Utilities/checkGeo.m",
        "Utilities/staticDetectorGeo.m",
        "Utilities/computeCOR.m",
        "Utilities/flatten_detector.m",
        "Utilities/plotgeometry.m",
    ],
    "utilities": [
        "Utilities/addCTnoise.m",
        "Utilities/im3DDenoise.m",
        "Utilities/im2DDenoise.m",
        "Utilities/im3Dnorm.m",
        "Utilities/plotImg.m",
        "Utilities/plotProj.m",
        "Utilities/plotSinogram.m",
        "Utilities/sheppLogan3D.m",
        "Utilities/cropCBCT.m",
        "Utilities/out_of_FOV_mask.m",
        "Utilities/Quality_measures/Measure_Quality.m",
        "Utilities/GPU/GpuIds.m",
        "Utilities/GPU/getGpuIds.m",
        "Utilities/GPU/getGpuNames.m",
        "Utilities/citeme.m",
    ],
    "io": [
        "Utilities/IO/Nikon/NikonDataLoader.m",
        "Utilities/IO/Bruker/BrukerDataLoader.m",
        "Utilities/IO/YXLON/YXLONDataLoader.m",
        "Utilities/IO/Diondo/DiondoDataLoader.m",
        "Utilities/IO/VarianCBCT/VarianDataLoader.m",
        "Utilities/IO/DXChange/DXChangeDataLoader.m",
        "Utilities/IO/Dicom/dicomDataLoader.m",
    ],
}

_LICENSE = re.compile(r"This file is part of the TIGRE|Copyright \(c\)|^-{5,}$")


def _files(entry):
    path = MATLAB / entry
    if path.is_dir():
        return sorted(p for p in path.glob("*.m"))
    return [path] if path.exists() else []


def parse(path):
    """Return (name, signature, help text) of a MATLAB function or class file."""
    lines = path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").split("\n")
    signature, i = "", 0
    while i < len(lines) and not lines[i].strip().startswith(("function", "classdef")):
        i += 1
    if i < len(lines):
        signature = re.sub(r"\s*%.*$", "", lines[i]).strip()
        i += 1
    # the help text: the first block of comment lines after the signature
    while i < len(lines) and not lines[i].strip():
        i += 1
    help_lines = []
    while i < len(lines) and lines[i].lstrip().startswith("%"):
        text = re.sub(r"^\s*%+ ?", "", lines[i]).rstrip()
        if _LICENSE.search(text):
            break
        help_lines.append(text)
        i += 1
    while help_lines and (not help_lines[-1].strip() or set(help_lines[-1].strip()) <= {"-", "%"}):
        help_lines.pop()
    return path.stem, signature, "\n".join(help_lines)


def _fence(text):
    fence = "```"
    while fence in text:
        fence += "`"
    return fence


def render(group):
    out = []
    for entry in GROUPS[group]:
        for path in _files(entry):
            name, signature, help_text = parse(path)
            rel = path.relative_to(MATLAB).as_posix()
            out.append(f"(matlab-{name.lower()})=\n## {name}\n")
            if signature:
                out.append(f"```matlab\n{signature}\n```\n")
            if help_text.strip():
                fence = _fence(help_text)
                out.append(f"{fence}text\n{help_text}\n{fence}\n")
            else:
                out.append("No help text.\n")
            out.append(f"Source: [`MATLAB/{rel}`]({SOURCE_URL}{rel})\n")
    return "\n".join(out)


def generate(app=None):
    OUT.mkdir(parents=True, exist_ok=True)
    for group in GROUPS:
        (OUT / f"{group}.md").write_text(render(group), encoding="utf-8", newline="\n")


def setup(app):
    app.connect("builder-inited", generate)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
