"""Render the gallery on a GPU machine and store it on the docs-gallery branch.

The website build on GitHub has no GPU, so it shows the gallery demos as code
only. This script runs the demos (through sphinx-gallery) on a computer with
TIGRE installed and a GPU, then saves the rendered gallery (pages, images,
animations) as a commit on the ``docs-gallery`` branch. The Docs workflow copies
that branch into ``docs/gallery`` before building, and sphinx-gallery reuses
each rendered page for as long as its demo is unchanged. A demo that changed
later appears as code only until the script is run again.

Run it from the root of an up-to-date checkout of master, in the Python
environment where TIGRE works:

    pip install -r docs/requirements.txt
    python docs/tools/build_gallery.py            # render, commit to docs-gallery
    git push --force origin docs-gallery          # publish
    # then: Actions > Docs > Run workflow (on master) to redeploy the site

Options:
    --fresh   re-run every demo instead of reusing the current docs-gallery pages
    --push    also push docs-gallery and start the Docs workflow (needs gh)

The script never touches your checked-out branch or working tree: it only
writes docs/gallery and docs/_build (both ignored by git) and moves the
docs-gallery branch.
"""
from __future__ import annotations

import argparse
import io
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
GALLERY = DOCS / "gallery"
BRANCH = "docs-gallery"


def git(*args: str, env: dict | None = None, capture: bool = True) -> str:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, env=env, check=True,
        stdout=subprocess.PIPE if capture else None, text=False,
    )
    return result.stdout.decode().strip() if capture and result.stdout else ""


def remote_branch_exists() -> bool:
    try:
        return bool(git("ls-remote", "--heads", "origin", BRANCH))
    except subprocess.CalledProcessError:
        return False


def restore_previous() -> None:
    """Unpack the current docs-gallery branch into docs/gallery."""
    if not remote_branch_exists():
        print(f"No {BRANCH} branch on origin yet: every demo will run.")
        return
    git("fetch", "origin", BRANCH, capture=False)
    archive = subprocess.run(
        ["git", "archive", "--format=zip", "FETCH_HEAD"], cwd=ROOT, check=True,
        stdout=subprocess.PIPE,
    ).stdout
    with zipfile.ZipFile(io.BytesIO(archive)) as zf:
        zf.extractall(GALLERY)
    print(f"Reusing the pages on origin/{BRANCH}; only changed demos will run.")


def check_tigre() -> None:
    try:
        import tigre  # noqa: F401
        from tigre.utilities import gpu
    except ImportError as err:
        sys.exit(f"TIGRE is not importable in this Python ({err}). Activate the right environment.")
    names = gpu.getGpuNames()
    if not names:
        sys.exit("TIGRE found no GPU.")
    print("GPUs:", ", ".join(names))


def build() -> None:
    env = dict(os.environ, TIGRE_GALLERY_RUN="1")
    cmd = [sys.executable, "-m", "sphinx", "-b", "html", str(DOCS), str(DOCS / "_build" / "html")]
    print(" ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)


def commit() -> str:
    """Commit docs/gallery as the new tip of the docs-gallery branch."""
    (GALLERY / ".gitattributes").write_text("* -text\n", encoding="utf-8")
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ, GIT_INDEX_FILE=str(pathlib.Path(tmp) / "index"))
        git("--work-tree", str(GALLERY), "add", "--all", "--force", ".", env=env)
        tree = git("write-tree", env=env)
    # One commit with no history: each render replaces the last, so old images
    # do not pile up in the repository. Publishing therefore needs --force.
    head = git("rev-parse", "--short", "HEAD")
    sha = git("commit-tree", tree, "-m", f"Rendered gallery from {head}")
    git("update-ref", f"refs/heads/{BRANCH}", sha)
    return sha


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--fresh", action="store_true", help="re-run every demo")
    parser.add_argument("--push", action="store_true", help="push and redeploy (needs gh)")
    args = parser.parse_args()

    check_tigre()
    if GALLERY.exists():
        shutil.rmtree(GALLERY)
    if not args.fresh:
        restore_previous()
    build()
    sha = commit()
    print(f"\n{BRANCH} is now {sha[:10]}.")
    if args.push:
        git("push", "--force", "origin", BRANCH, capture=False)
        subprocess.run(["gh", "workflow", "run", "docs.yml", "--ref", "master"], cwd=ROOT, check=False)
    else:
        print(f"Publish with:  git push --force origin {BRANCH}")
        print("Then redeploy: Actions > Docs > Run workflow (branch master).")


if __name__ == "__main__":
    main()
