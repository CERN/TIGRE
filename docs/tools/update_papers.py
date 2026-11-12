"""Fill in metadata for docs/papers.yaml from Semantic Scholar.

Entries in docs/papers.yaml can be as short as an identifier and an area:

    - doi: "10.1088/1361-6560/ae8123"
      area: medical

This script looks up every entry that is missing a title, authors, venue or
year, and fills those fields in (full author lists). Fields already present are
never overwritten unless --refresh is given, so manual corrections stick.

Usage:
    python docs/tools/update_papers.py               # fill in missing metadata
    python docs/tools/update_papers.py --refresh     # re-fetch everything
    python docs/tools/update_papers.py --candidates  # list papers citing TIGRE
                                                     # that are not in the file yet

Identifiers: ``doi``, ``arxiv`` (e.g. "2405.20693") or ``s2`` (a Semantic
Scholar paper id). Set S2_API_KEY to use a Semantic Scholar API key.

Network access is needed, so this runs in the "Papers list" GitHub workflow or
on a maintainer's machine, never during the docs build.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
import urllib.error
import urllib.request

import yaml

PAPERS = pathlib.Path(__file__).resolve().parent.parent / "papers.yaml"
API = "https://api.semanticscholar.org/graph/v1"
FIELDS = "title,authors,year,venue,journal,externalIds,publicationVenue"
TIGRE_PAPERS = [
    "DOI:10.1088/2057-1976/2/5/055010",  # 2016, BPEX
    "DOI:10.1016/j.jpdc.2020.07.004",  # 2020, JPDC
    "DOI:10.1088/2631-8695/adbb3a",  # 2025, TIGRE v3
]
KEY_ORDER = ["title", "authors", "venue", "year", "doi", "url", "arxiv", "s2",
             "area", "note", "verified"]
HEADER = """\
# Published work that USES TIGRE (not only cites it).
#
# This file generates docs/papers.md. To add a paper, append an entry with an
# identifier and an area, then open a pull request:
#
#   - doi: "10.xxxx/yyyy"        # or  arxiv: "2405.20693"
#     area: medical
#
# The "Papers list" workflow (or `python docs/tools/update_papers.py`) fills in
# title, authors, venue and year. Fields you write yourself are never overwritten.
#
# Fields: title, authors, venue, year | doi, url, arxiv, s2 | area | note
#   area: medical, particle, industrial, synchrotron, methods or other
#   note: optional, one line on how TIGRE was used
#   verified: how we know the paper uses TIGRE (e.g. citation-context: the
#     sentence citing TIGRE says it was used; author: confirmed by the authors)
# Papers about TIGRE itself are listed on the Citing page, not here.
"""


def _request(url, body=None):
    headers = {"User-Agent": "TIGRE-docs (github.com/CERN/TIGRE)"}
    if os.environ.get("S2_API_KEY"):
        headers["x-api-key"] = os.environ["S2_API_KEY"]
    data = None
    if body is not None:
        data = json.dumps(body).encode()
        headers["Content-Type"] = "application/json"
    for attempt in range(6):
        req = urllib.request.Request(url, data=data, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as err:
            if err.code in (429, 500, 502, 503, 504) and attempt < 5:
                time.sleep(5 * (attempt + 1))
                continue
            raise


def _s2_id(entry):
    if entry.get("s2"):
        return entry["s2"]
    if entry.get("doi"):
        return f"DOI:{entry['doi']}"
    if entry.get("arxiv"):
        return f"ARXIV:{entry['arxiv']}"
    return None


def _venue(paper):
    journal = paper.get("journal") or {}
    name = journal.get("name") or (paper.get("publicationVenue") or {}).get("name") or paper.get("venue")
    if name in (None, "", "arXiv.org", "ArXiv"):
        arxiv = (paper.get("externalIds") or {}).get("ArXiv")
        return f"arXiv:{arxiv}" if arxiv else "Preprint"
    extra = ", ".join(x.strip() for x in (journal.get("volume"), journal.get("pages")) if x and x.strip())
    return f"{name} {extra}".strip() if extra else name


def _fill(entry, paper, refresh):
    ext = paper.get("externalIds") or {}
    new = {
        "title": (paper.get("title") or "").strip() or None,
        "authors": [a["name"] for a in paper.get("authors") or []] or None,
        "venue": _venue(paper),
        "year": paper.get("year"),
    }
    for key, value in new.items():
        if value is not None and (refresh or not entry.get(key)):
            entry[key] = value
    doi = ext.get("DOI")
    if doi and not doi.startswith("10.48550") and not entry.get("doi"):
        entry["doi"] = doi
    if not entry.get("doi") and not entry.get("url") and ext.get("ArXiv"):
        entry["url"] = f"https://arxiv.org/abs/{ext['ArXiv']}"


def _ordered(entry):
    keys = [k for k in KEY_ORDER if k in entry] + [k for k in entry if k not in KEY_ORDER]
    return {k: entry[k] for k in keys}


def _dump(entries):
    blocks = []
    for e in entries:
        text = yaml.safe_dump([_ordered(e)], allow_unicode=True, sort_keys=False, width=4096, default_flow_style=None)
        blocks.append(text)
    return HEADER + "\n" + "\n".join(blocks)


def update(refresh=False):
    entries = yaml.safe_load(PAPERS.read_text(encoding="utf-8")) or []
    todo = [e for e in entries if _s2_id(e) and (refresh or not all(e.get(k) for k in ("title", "authors", "venue", "year")))]
    print(f"{len(entries)} entries, {len(todo)} to look up")
    for start in range(0, len(todo), 400):
        chunk = todo[start:start + 400]
        papers = _request(f"{API}/paper/batch?fields={FIELDS}", {"ids": [_s2_id(e) for e in chunk]})
        for entry, paper in zip(chunk, papers):
            if paper is None:
                print(f"  not found: {_s2_id(entry)}", file=sys.stderr)
                continue
            _fill(entry, paper, refresh)
    PAPERS.write_text(_dump(entries), encoding="utf-8")
    missing = [_s2_id(e) or e.get("title") for e in entries if not all(e.get(k) for k in ("title", "authors", "venue", "year"))]
    if missing:
        print("Still incomplete (fix by hand):", *missing, sep="\n  ", file=sys.stderr)
        return 1
    return 0


def candidates():
    """Print papers citing TIGRE that are not in papers.yaml, with the citing sentences."""
    entries = yaml.safe_load(PAPERS.read_text(encoding="utf-8")) or []
    known = {str(e.get(k)).lower() for e in entries for k in ("doi", "arxiv", "s2", "title") if e.get(k)}
    seen = {}
    for tigre in TIGRE_PAPERS:
        offset = 0
        while offset is not None:
            page = _request(f"{API}/paper/{tigre}/citations?fields=title,year,externalIds,contexts&limit=100&offset={offset}")
            for item in page["data"]:
                paper = item["citingPaper"]
                if paper.get("paperId"):
                    seen.setdefault(paper["paperId"], (paper, []))[1].extend(item.get("contexts") or [])
            offset = page.get("next")
            time.sleep(1)
    print("| Year | Paper | Sentences citing TIGRE |\n|---|---|---|")
    for pid, (paper, contexts) in sorted(seen.items(), key=lambda kv: -(kv[1][0].get("year") or 0)):
        ext = paper.get("externalIds") or {}
        ids = {pid.lower(), (paper.get("title") or "").lower(), str(ext.get("DOI", "")).lower(), str(ext.get("ArXiv", "")).lower()}
        if ids & known:
            continue
        ctx = " / ".join(c.replace("|", " ") for c in contexts) or "(no citation context available)"
        print(f"| {paper.get('year')} | {paper.get('title')} (s2: {pid}) | {ctx} |")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--refresh", action="store_true", help="re-fetch metadata for every entry")
    parser.add_argument("--candidates", action="store_true", help="list citing papers not yet in the file")
    args = parser.parse_args()
    sys.exit(candidates() if args.candidates else update(args.refresh))
