"""Generate the "Papers using TIGRE" list from docs/papers.yaml.

Writes ``docs/_generated/papers.md``, grouped by application area and sorted
newest first. ``papers.md`` includes it. An unknown area fails the build.
Entries without title/authors/venue/year yet (added with only an identifier)
are skipped until ``docs/tools/update_papers.py`` fills them in.
"""

from __future__ import annotations

import pathlib

import yaml
from sphinx.errors import ExtensionError
from sphinx.util import logging

logger = logging.getLogger(__name__)

DOCS = pathlib.Path(__file__).resolve().parent.parent
SRC = DOCS / "papers.yaml"
OUT = DOCS / "_generated" / "papers.md"

AREAS = {
    "medical": "Medical imaging",
    "particle": "Proton and ion CT",
    "industrial": "Industrial and non-destructive testing",
    "synchrotron": "Synchrotron and neutron imaging",
    "materials": "X-ray physics, detectors and materials",
    "methods": "Reconstruction methods and machine learning",
    "other": "Other applications",
}
REQUIRED = ("title", "authors", "venue", "year")


def _authors(names):
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def _entry(p):
    title = p["title"].rstrip(".")
    if p.get("doi"):
        link = f"https://doi.org/{p['doi']}"
    elif p.get("url"):
        link = p["url"]
    elif p.get("arxiv"):
        link = f"https://arxiv.org/abs/{p['arxiv']}"
    else:
        link = None
    title_md = f"[{title}]({link})" if link else title
    line = f"- {_authors(p['authors'])}. **{title_md}.** *{p['venue']}*, {p['year']}."
    if p.get("note"):
        line += f"\n  {p['note']}"
    return line


def render(papers: list[dict]) -> str:
    complete = []
    for i, p in enumerate(papers):
        if p.get("area") not in AREAS:
            raise ExtensionError(
                f"papers.yaml entry {i + 1}: unknown area {p.get('area')!r}; "
                f"use one of {', '.join(AREAS)}"
            )
        missing = [k for k in REQUIRED if not p.get(k)]
        if missing:
            # Entries added with only an identifier are filled in by
            # docs/tools/update_papers.py; until then they are left out.
            ident = p.get("doi") or p.get("arxiv") or p.get("s2") or p.get("title")
            logger.info(f"papers.yaml: skipping {ident} (missing {', '.join(missing)})")
            continue
        complete.append(p)
    papers = complete
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
