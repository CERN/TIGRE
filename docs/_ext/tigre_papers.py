"""Generate the "Papers using TIGRE" list from docs/papers.yaml.

Writes ``docs/_generated/papers.md``, grouped by application area and sorted
newest first. ``papers.md`` includes it. Unknown areas or missing required
fields fail the build, so mistakes in a pull request are caught by CI.
"""

from __future__ import annotations

import pathlib

import yaml
from sphinx.errors import ExtensionError

DOCS = pathlib.Path(__file__).resolve().parent.parent
SRC = DOCS / "papers.yaml"
OUT = DOCS / "_generated" / "papers.md"

AREAS = {
    "medical": "Medical imaging",
    "particle": "Proton and ion CT",
    "industrial": "Industrial and non-destructive testing",
    "synchrotron": "Synchrotron and neutron imaging",
    "methods": "Reconstruction methods and machine learning",
    "other": "Other applications",
}
REQUIRED = ("title", "authors", "venue", "year", "area")


def _authors(names):
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def _entry(p):
    title = p["title"].rstrip(".")
    link = f"https://doi.org/{p['doi']}" if p.get("doi") else p.get("url")
    title_md = f"[{title}]({link})" if link else title
    line = f"- {_authors(p['authors'])}. **{title_md}.** *{p['venue']}*, {p['year']}."
    if p.get("note"):
        line += f"\n  {p['note']}"
    return line


def render(papers: list[dict]) -> str:
    for i, p in enumerate(papers):
        missing = [k for k in REQUIRED if k not in p]
        if missing:
            raise ExtensionError(f"papers.yaml entry {i + 1}: missing {', '.join(missing)}")
        if p["area"] not in AREAS:
            raise ExtensionError(
                f"papers.yaml entry {i + 1}: unknown area {p['area']!r}; "
                f"use one of {', '.join(AREAS)}"
            )
    parts = ["<!-- Generated from docs/papers.yaml by docs/_ext/tigre_papers.py. Do not edit. -->", ""]
    for key, heading in AREAS.items():
        group = sorted((p for p in papers if p["area"] == key), key=lambda p: -int(p["year"]))
        if not group:
            continue
        parts += [f"## {heading}", ""] + [_entry(p) for p in group] + [""]
    return "\n".join(parts)


def _generate(app):
    papers = yaml.safe_load(SRC.read_text(encoding="utf-8")) or []
    OUT.parent.mkdir(exist_ok=True)
    text = render(papers)
    if not OUT.exists() or OUT.read_text(encoding="utf-8") != text:
        OUT.write_text(text, encoding="utf-8")


def setup(app):
    app.connect("builder-inited", _generate)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
