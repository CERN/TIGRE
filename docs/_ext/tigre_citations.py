"""Generate the citation list on the Citing page from CITATION.cff.

CITATION.cff is the single source of truth for how to cite TIGRE. At build
time this extension writes ``docs/_generated/citations.md`` (a MyST fragment)
with a formatted reference and a BibTeX entry for ``preferred-citation`` and
for every entry under ``references``, in that order. ``citing.md`` includes that file.
"""

from __future__ import annotations

import pathlib
import re
import unicodedata

import yaml

DOCS = pathlib.Path(__file__).resolve().parent.parent
ROOT = DOCS.parent
CFF = ROOT / "CITATION.cff"
OUT = DOCS / "_generated" / "citations.md"


def _names(authors):
    return [a.get("name") or f"{a['given-names']} {a['family-names']}" for a in authors]


def _initials(given):
    parts = re.split(r"([\s-])", given)
    return "".join(p[0].upper() + "." if p.strip(" -") else p for p in parts if p)


def _short_names(authors):
    out = []
    for a in authors:
        if "family-names" in a:
            out.append(f"{_initials(a['given-names'])} {a['family-names']}")
        else:
            out.append(a["name"])
    return out


def _join(names):
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def _ascii(text):
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()


def _bibkey(ref):
    first = ref["authors"][0].get("family-names", "tigre")
    word = re.sub(r"[^A-Za-z]", "", ref["title"].split()[0])
    return f"{_ascii(first)}{ref['year']}{word}".replace(" ", "")


def _formatted(ref):
    text = f"{_join(_short_names(ref['authors']))}. *{ref['title']}.* {ref['journal']}"
    if "volume" in ref:
        text += f" **{ref['volume']}**"
        if "issue" in ref:
            text += f"({ref['issue']})"
    if "start" in ref:
        text += f", {ref['start']}"
    text += f" ({ref['year']})."
    if "doi" in ref:
        text += f" [doi:{ref['doi']}](https://doi.org/{ref['doi']})"
    return text


def _bibtex(ref):
    fields = [
        ("author", " and ".join(
            f"{a['family-names']}, {a['given-names']}" if "family-names" in a else a["name"]
            for a in ref["authors"]
        )),
        ("title", "{" + ref["title"] + "}"),
        ("journal", ref.get("journal")),
        ("volume", ref.get("volume")),
        ("number", ref.get("issue")),
        ("pages", ref.get("start")),
        ("year", ref.get("year")),
        ("doi", ref.get("doi")),
    ]
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields if v is not None)
    return f"@article{{{_bibkey(ref)},\n{body}\n}}"


def _block(ref, heading):
    return "\n".join([
        f"## {heading}",
        "",
        _formatted(ref),
        "",
        ":::{dropdown} BibTeX",
        "```bibtex",
        _bibtex(ref),
        "```",
        ":::",
        "",
    ])


def _heading(ref):
    first = ref["authors"][0].get("family-names", ref["authors"][0].get("name"))
    etal = " et al." if len(ref["authors"]) > 1 else ""
    return f"{first}{etal} {ref['year']}, {ref['journal']}"


def render(cff: dict) -> str:
    parts = [
        "<!-- Generated from CITATION.cff by docs/_ext/tigre_citations.py. Do not edit. -->",
        "",
        _block(cff["preferred-citation"], _heading(cff["preferred-citation"])),
    ]
    for ref in cff.get("references", []):
        parts.append(_block(ref, _heading(ref)))
    return "\n".join(parts)


def _generate(app):
    cff = yaml.safe_load(CFF.read_text(encoding="utf-8"))
    OUT.parent.mkdir(exist_ok=True)
    text = render(cff)
    if not OUT.exists() or OUT.read_text(encoding="utf-8") != text:
        OUT.write_text(text, encoding="utf-8")


def setup(app):
    app.connect("builder-inited", _generate)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
