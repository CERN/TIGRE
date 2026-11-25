"""MATLAB API reference pages, from the help text of TIGRE's .m files.

For each group in ``GROUPS`` this writes ``docs/_generated/matlab/<group>.md``
(a MyST fragment), which ``docs/api/matlab/<group>.md`` includes. Each function
gets a heading, its signature, its help text (what ``help NAME`` prints in
MATLAB, without the license block) laid out as text and an options list, the
original help text in a drop-down, and a link to the source. A table at the
top of each page lists the functions with their first sentence.
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
        "Utilities/tigreVersion.m",
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


_OPTION = re.compile(r"^(\s*)(?:\*\s*)?'([A-Za-z_][\w]*)'\s*:?\s*(.*)$")


def _escape(text):
    """Make plain help text safe as Markdown."""
    text = text.replace("\\", "\\\\").replace("*", "\\*").replace("<", "\\<")
    text = text.replace("|", "\\|").replace("`", "\\`")
    return text


def _code_calls(text, name):
    """Put calls of the function itself, like OS_SART(PROJ,GEO), in code."""
    return re.sub(rf"\b({re.escape(name)}\([^()]*\))", r"`\1`", text, flags=re.I)


def _paragraphs(lines):
    paras, cur = [], []
    for line in lines:
        if line.strip():
            cur.append(line.strip())
        elif cur:
            paras.append(" ".join(cur))
            cur = []
    if cur:
        paras.append(" ".join(cur))
    return paras


def split_help(name, help_text):
    """Split MATLAB help text into intro paragraphs and (option, description) pairs."""
    lines = help_text.split("\n")
    matches = [(k, _OPTION.match(l)) for k, l in enumerate(lines)]
    matches = [(k, m) for k, m in matches if m]
    if not matches:
        return _paragraphs(lines), []
    # The first option sets the indent; deeper quoted words are values of an option.
    indent = len(matches[0][1].group(1))
    starts = [k for k, m in matches if len(m.group(1)) <= indent + 2]
    intro = _paragraphs(lines[: starts[0]])
    options = []
    for n, k in enumerate(starts):
        end = starts[n + 1] if n + 1 < len(starts) else len(lines)
        m = _OPTION.match(lines[k])
        body = [m.group(3)] + lines[k + 1 : end]
        options.append((m.group(2), body))
    return intro, options


def _option_description(body):
    """The description of one option: a paragraph, then any sub-values as a list."""
    text, items = [], []
    for line in body:
        sub = _OPTION.match(line) if line.strip() else None
        if sub and (text or items):
            items.append([sub.group(2), sub.group(3).strip()])
        elif items and line.strip():
            items[-1][1] += " " + line.strip()
        elif line.strip():
            text.append(line.strip())
    out = _escape(" ".join(text))
    for value, desc in items:
        desc = desc.lstrip(":").strip()
        out += f"\n  - `'{value}'`: {_escape(desc)}"
    return out


def summary(name, help_text):
    """One sentence about the function, for the table at the top of the page."""
    intro, _ = split_help(name, help_text)
    if not intro:
        return ""
    first = intro[0]
    first = re.sub(rf"^(geo\s*=\s*)?{re.escape(name)}(\([^()]*\))?\s*", "", first, flags=re.I)
    first = re.split(r"(?<=\.)\s", first)[0].strip()
    if len(first) > 140:
        first = first[:137].rsplit(" ", 1)[0] + "..."
    return _escape(first[:1].upper() + first[1:]) if first else ""


def render(group):
    entries = []
    for entry in GROUPS[group]:
        for path in _files(entry):
            entries.append((path, *parse(path)))
    out = ["| Function | What it does |", "|---|---|"]
    for _, name, _, help_text in entries:
        out.append(f"| {{ref}}`{name} <matlab-{name.lower()}>` | {summary(name, help_text)} |")
    out.append("")
    for path, name, signature, help_text in entries:
        rel = path.relative_to(MATLAB).as_posix()
        out.append(f"(matlab-{name.lower()})=\n## {name}\n")
        if signature:
            out.append(f"```matlab\n{signature}\n```\n")
        if not help_text.strip():
            out.append("No help text.\n")
        else:
            intro, options = split_help(name, help_text)
            for para in intro:
                out.append(_code_calls(_escape(para), name) + "\n")
            if options:
                if not (intro and intro[-1].rstrip().endswith(":")):
                    out.append("Options:\n")
                for option, body in options:
                    out.append(f"`'{option}'`\n: {_option_description(body)}\n")
            fence = _fence(help_text)
            out.append(f":::{{dropdown}} Help text as in MATLAB\n{fence}text\n{help_text}\n{fence}\n:::\n")
        out.append(f"Source: [`MATLAB/{rel}`]({SOURCE_URL}{rel})\n")
    return "\n".join(out)


def generate(app=None):
    OUT.mkdir(parents=True, exist_ok=True)
    for group in GROUPS:
        (OUT / f"{group}.md").write_text(render(group), encoding="utf-8", newline="\n")


def setup(app):
    app.connect("builder-inited", generate)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
