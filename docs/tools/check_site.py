"""Weekly checks of the website: broken links and out-of-date demo renders.

Run by .github/workflows/site-checks.yml after a linkcheck build, with the
docs-demos branch unpacked into docs/demos. Writes a Markdown report and
prints the number of problems found (0 when everything is fine):

    python docs/tools/check_site.py --report report.md

Checks:

- **Broken links**: the "broken" entries of Sphinx's linkcheck output
  (docs/_build/linkcheck/output.json). Redirects and timeouts are not counted.
- **Stale demos**: a demo whose source changed since it was last rendered with
  docs/tools/build_demos.py. Its page then shows the code only, without
  outputs, until the demos are rendered again. This compares each demo's
  md5 with the one stored on the docs-demos branch, as sphinx-gallery does.
- **Old Colab notebooks**: notebooks whose first cell installs an older
  TIGRE than the version in pyproject.toml (after a release, before the demos
  are rendered again with ``build_demos.py --fresh``).
"""
import argparse
import hashlib
import json
import pathlib
import re
import sys
import tempfile

DOCS = pathlib.Path(__file__).resolve().parents[1]
ROOT = DOCS.parent
SITE = "https://cern.github.io/TIGRE/"


def broken_links(output_json: pathlib.Path) -> list[dict]:
    if not output_json.exists():
        return []
    broken = []
    for line in output_json.read_text(encoding="utf-8").splitlines():
        entry = json.loads(line)
        if entry.get("status") == "broken":
            broken.append(entry)
    return broken


def _demo_sources() -> tuple[pathlib.Path, set[str]]:
    sys.path.insert(0, str(DOCS / "_ext"))
    import tigre_demos  # noqa: E402

    out = tigre_demos.generate(pathlib.Path(tempfile.mkdtemp()) / "demos_src")
    return out, set(tigre_demos.NOT_RUN)


def stale_demos(rendered: pathlib.Path) -> list[str]:
    """Rendered demos whose source changed (or which were never rendered)."""
    sources, not_run = _demo_sources()
    stale = []
    for src in sorted(sources.rglob("*.py")):
        if src.name in not_run:
            continue
        rel = src.relative_to(sources)
        md5 = rendered / (str(rel) + ".md5")
        current = hashlib.md5(src.read_bytes()).hexdigest()
        if not md5.exists() or md5.read_text(encoding="utf-8").strip() != current:
            stale.append(rel.with_suffix("").as_posix())
    return stale


def old_notebooks(rendered: pathlib.Path) -> tuple[str, list[str]]:
    """Notebooks whose Colab install cell is not for the current version."""
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    version = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M).group(1)
    old = []
    for nb in sorted(rendered.rglob("*.ipynb")):
        cells = json.loads(nb.read_text(encoding="utf-8")).get("cells", [])
        first = "".join(cells[0].get("source", "")) if cells else ""
        if f"releases/download/v{version}/" not in first:
            old.append(nb.relative_to(rendered).with_suffix("").as_posix())
    return version, old


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--linkcheck", type=pathlib.Path,
                        default=DOCS / "_build" / "linkcheck" / "output.json")
    parser.add_argument("--demos", type=pathlib.Path, default=DOCS / "demos",
                        help="the docs-demos branch, unpacked")
    args = parser.parse_args()

    links = broken_links(args.linkcheck)
    have_renders = args.demos.is_dir() and any(args.demos.rglob("*.md5"))
    stale = stale_demos(args.demos) if have_renders else []
    version, old = old_notebooks(args.demos) if have_renders else ("", [])

    lines = [
        "This issue is opened and updated by the weekly "
        "[site checks](../blob/master/.github/workflows/site-checks.yml) "
        "and closed when they pass.",
        "",
    ]
    if links:
        lines += [f"## Broken links ({len(links)})", "",
                  "| Page | Link | Error |", "|---|---|---|"]
        for e in links:
            page = SITE + e["filename"].replace(".md", ".html").replace(".rst", ".html")
            info = str(e.get("info", "")).replace("|", "\\|").replace("\n", " ")[:150]
            lines.append(f"| [{e['filename']}]({page}) (line {e.get('lineno', '?')}) "
                         f"| {e['uri']} | {info} |")
        lines += ["", "Fix or remove each link, or add it to `linkcheck_ignore` in "
                  "`docs/conf.py` if the site only blocks automated checkers.", ""]
    if stale:
        lines += [f"## Demos to render again ({len(stale)})", "",
                  "Their code changed since they were last rendered, so their pages show "
                  "the code without outputs:", ""]
        lines += [f"- `{d}`" for d in stale]
        lines += ["", "Run `python docs/tools/build_demos.py` on a GPU machine and push "
                  "`docs-demos` (see the script's docstring).", ""]
    if old:
        lines += [f"## Colab notebooks for an older TIGRE ({len(old)})", "",
                  f"Their first cell does not install TIGRE {version}:", ""]
        lines += [f"- `{d}`" for d in old[:10]]
        if len(old) > 10:
            lines.append(f"- and {len(old) - 10} more")
        lines += ["", "Run `python docs/tools/build_demos.py --fresh` on a GPU machine "
                  "and push `docs-demos`.", ""]
    if not have_renders:
        lines += ["The docs-demos branch was not found, so the demos were not checked.", ""]

    problems = len(links) + len(stale) + len(old)
    args.report.write_text("\n".join(lines), encoding="utf-8")
    print(problems)
    return 0


if __name__ == "__main__":
    sys.exit(main())
