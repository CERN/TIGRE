"""Build the review queue for "Papers using TIGRE".

Collects every paper citing the TIGRE papers (Semantic Scholar) that is not yet
in docs/papers.yaml or docs/tools/papers_reviewed.yaml, and gathers the
sentences that mention TIGRE:

1. citation contexts from Semantic Scholar, where available;
2. otherwise, sentences extracted from the open-access full text (PDF).

Papers that OpenAlex lists as citing TIGRE but Semantic Scholar does not are
added too (OpenAlex has no citation contexts, so those rely on full text).

Writes docs/tools/review_queue.json for the review page. Runs in CI (needs
network): see .github/workflows/papers-review.yml.
"""

from __future__ import annotations

import io
import json
import pathlib
import re
import sys
import time
import urllib.request

import yaml

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from update_papers import API, PAPERS, TIGRE_PAPERS, _request  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
REVIEWED = HERE / "papers_reviewed.yaml"
HINTS = HERE / "review_hints.json"
OUT = HERE / "review_queue.json"
MENTION = re.compile(r"TIGRE|TIGER\b|TI-GRE|Tomographic Iterative GPU", re.I)
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/130 Safari/537.36"


def guess_area(title, venue, summary):
    s = " ".join(x for x in (title, venue, summary) if x)
    rules = [
        ("particle", r"proton|\bion (ct|imaging|computed)|carbon ion|hadron"),
        ("synchrotron", r"synchrotron|ptycho|neutron|beamline"),
        ("industrial", r"laminograph|batter|industrial|non-?destructive|nondestructive|metrolog|additive manufactur|composite|weld|casting|precision eng"),
        ("materials", r"photon-counting detector|detector calibration|spectral (x-ray|imaging)|material decomposition|scintillat|grating interferom|hyperspectral|walnut|rock|geolog|soil|wood"),
        ("methods", r"gaussian splat|nerf|neural|diffusion|implicit|deep|learning|network|transformer|algorithm|regulari[sz]|optimi[sz]ation|iterative reconstruction|inverse problem|reconstruction method|framework"),
        ("medical", r"patient|clinical|radiother|dental|breast|lung|cardiac|coronary|brain|spine|bone|tumou?r|dose|mammograph|cbct|angiograph|surg"),
    ]
    for area, pattern in rules:
        if re.search(pattern, s, re.I):
            return area
    return "other"


def sentences(text):
    text = re.sub(r"-\n(?=[a-z])", "", text)
    text = re.sub(r"\s+", " ", text)
    return re.split(r"(?<=[.!?])\s+(?=[A-Z(\[])", text)


def fulltext_sentences(url):
    try:
        from pypdf import PdfReader
    except ImportError:
        return None, "pypdf missing"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/pdf,*/*"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = resp.read(40_000_000)
        if not data.startswith(b"%PDF"):
            return None, "not a PDF (probably a landing page)"
        reader = PdfReader(io.BytesIO(data))
        text = "\n".join((page.extract_text() or "") for page in reader.pages)
    except Exception as err:  # network errors, paywalls, broken PDFs
        return None, f"{type(err).__name__}"
    # Drop the reference list: it always mentions TIGRE and says nothing about use.
    cut = re.search(r"\n\s*(References|REFERENCES|Bibliography)\s*\n", text)
    body = text[: cut.start()] if cut else text
    found = [s.strip() for s in sentences(body) if MENTION.search(s) and 20 < len(s) < 700]
    return found[:8], "ok" if found else "no mention in body text"


OPENALEX = "https://api.openalex.org"
MAILTO = "tigre.toolbox@gmail.com"


def _get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": f"TIGRE-docs (mailto:{MAILTO})"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except Exception:
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))


def _abstract(inv):
    if not inv:
        return ""
    words = sorted((pos, w) for w, positions in inv.items() for pos in positions)
    return " ".join(w for _, w in words)


def openalex_citing():
    """Works citing the TIGRE papers according to OpenAlex."""
    ids = []
    for tigre in TIGRE_PAPERS:
        work = _get_json(f"{OPENALEX}/works/doi:{tigre.removeprefix('DOI:')}?mailto={MAILTO}")
        ids.append(work["id"].rsplit("/", 1)[-1])
    out, cursor = [], "*"
    while cursor:
        page = _get_json(
            f"{OPENALEX}/works?filter=cites:{'|'.join(ids)}&per-page=200&cursor={cursor}&mailto={MAILTO}"
            "&select=id,doi,title,publication_year,primary_location,best_oa_location,abstract_inverted_index,ids"
        )
        out.extend(page["results"])
        cursor = page["meta"].get("next_cursor")
        time.sleep(0.5)
    return out


def known_ids():
    ids = set()
    for path in (PAPERS, REVIEWED):
        if path.exists():
            for e in yaml.safe_load(path.read_text(encoding="utf-8")) or []:
                for k in ("doi", "arxiv", "s2"):
                    if e.get(k):
                        ids.add(str(e[k]).lower())
    return ids


def main():
    hints = json.loads(HINTS.read_text()) if HINTS.exists() else {}
    hint_of = {pid: label for label, pids in hints.items() for pid in pids}
    known = known_ids()
    fields = "title,year,venue,externalIds,openAccessPdf,abstract,tldr,contexts"
    seen = {}
    for tigre in TIGRE_PAPERS:
        offset = 0
        while offset is not None:
            page = _request(f"{API}/paper/{tigre}/citations?fields={fields}&limit=100&offset={offset}")
            for item in page["data"]:
                p = item["citingPaper"]
                if not p.get("paperId"):
                    continue
                entry = seen.setdefault(p["paperId"], {"paper": p, "contexts": [], "cites": []})
                entry["cites"].append(tigre)
                for c in item.get("contexts") or []:
                    if c not in entry["contexts"]:
                        entry["contexts"].append(c)
            offset = page.get("next")
            time.sleep(1.5)

    queue = []
    for pid, entry in seen.items():
        p = entry["paper"]
        ext = p.get("externalIds") or {}
        ids = {pid.lower(), str(ext.get("DOI", "")).lower(), str(ext.get("ArXiv", "")).lower()}
        if ids & known or hint_of.get(pid) == "own":
            continue
        doi = ext.get("DOI")
        arxiv = ext.get("ArXiv")
        link = f"https://doi.org/{doi}" if doi else (f"https://arxiv.org/abs/{arxiv}" if arxiv else f"https://www.semanticscholar.org/paper/{pid}")
        summary = ((p.get("tldr") or {}).get("text") or (p.get("abstract") or "")[:400]).strip()
        item = {
            "s2": pid, "doi": doi, "arxiv": arxiv, "title": p.get("title"), "year": p.get("year"),
            "venue": p.get("venue"), "link": link, "summary": summary,
            "area": guess_area(p.get("title"), p.get("venue"), summary),
            "claude": hint_of.get(pid), "sentences": entry["contexts"], "source": "citation-context",
            "origin": "semanticscholar",
        }
        if not item["sentences"]:
            pdf = (p.get("openAccessPdf") or {}).get("url") or (f"https://arxiv.org/pdf/{arxiv}" if arxiv else None)
            if pdf:
                found, status = fulltext_sentences(pdf)
                item["source"] = "full-text" if found else f"full-text: {status}"
                item["sentences"] = found or []
                item["pdf"] = pdf
                time.sleep(1)
            else:
                item["source"] = "no open-access text"
        queue.append(item)

    # Papers OpenAlex knows about but Semantic Scholar does not.
    s2_dois = {str((e["paper"].get("externalIds") or {}).get("DOI", "")).lower() for e in seen.values()}
    s2_titles = {re.sub(r"\W+", "", (e["paper"].get("title") or "").lower()) for e in seen.values()}
    try:
        oa_works = openalex_citing()
    except Exception as err:
        print("OpenAlex failed:", err)
        oa_works = []
    for w in oa_works:
        doi = (w.get("doi") or "").removeprefix("https://doi.org/") or None
        tkey = re.sub(r"\W+", "", (w.get("title") or "").lower())
        if (doi and doi.lower() in s2_dois) or tkey in s2_titles or (doi and doi.lower() in known):
            continue
        src = (w.get("primary_location") or {}).get("source") or {}
        summary = _abstract(w.get("abstract_inverted_index"))[:400]
        item = {
            "s2": None, "openalex": w["id"].rsplit("/", 1)[-1], "doi": doi, "arxiv": None,
            "title": w.get("title"), "year": w.get("publication_year"), "venue": src.get("display_name"),
            "link": w.get("doi") or w["id"], "summary": summary,
            "area": guess_area(w.get("title"), src.get("display_name"), summary),
            "claude": None, "sentences": [], "source": "no open-access text",
        }
        pdf = (w.get("best_oa_location") or {}).get("pdf_url")
        if pdf:
            found, status = fulltext_sentences(pdf)
            item["source"] = "full-text" if found else f"full-text: {status}"
            item["sentences"] = found or []
            item["pdf"] = pdf
            time.sleep(1)
        item["origin"] = "openalex"
        queue.append(item)
    print(len(oa_works), "OpenAlex citing works")

    order = {"unclear": 0, None: 1, "mention": 2}
    queue.sort(key=lambda q: (not q["sentences"], order.get(q["claude"], 1), -(q["year"] or 0)))
    OUT.write_text(json.dumps({"generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()), "papers": queue}, ensure_ascii=False, indent=1), encoding="utf-8")
    stats = {}
    for q in queue:
        key = q["source"] if q["sentences"] or not q["source"].startswith("full-text:") else "full-text failed"
        stats[key] = stats.get(key, 0) + 1
    print(len(queue), "papers in queue", stats)


if __name__ == "__main__":
    main()
