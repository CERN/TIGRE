"""Render the demos on a GPU machine and store it on the docs-demos branch.

The website build on GitHub has no GPU, so it shows the demos as code
only. This script runs the demos (through sphinx-gallery) on a computer with
TIGRE installed and a GPU, then saves the rendered demos (pages, images,
animations) as a commit on the ``docs-demos`` branch. The Docs workflow copies
that branch into ``docs/demos`` before building, and sphinx-gallery reuses
each rendered page for as long as its demo is unchanged. A demo that changed
later appears as code only until the script is run again.

Run it from the root of an up-to-date checkout of master, in the Python
environment where TIGRE works:

    pip install -r docs/requirements.txt
    python docs/tools/build_demos.py            # render, commit to docs-demos
    git push --force origin docs-demos          # publish
    # then: Actions > Docs > Run workflow (on master) to redeploy the site

Options:
    --fresh   re-run every demo instead of reusing the current docs-demos pages
    --push    also push docs-demos and start the Docs workflow (needs gh)

The script never touches your checked-out branch or working tree: it only
writes docs/demos and docs/_build (both ignored by git) and moves the
docs-demos branch.
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
DEMOS = DOCS / "demos"
BRANCH = "docs-demos"


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
    """Unpack the current docs-demos branch into docs/demos."""
    if not remote_branch_exists():
        print(f"No {BRANCH} branch on origin yet: every demo will run.")
        return
    git("fetch", "origin", BRANCH, capture=False)
    archive = subprocess.run(
        ["git", "archive", "--format=zip", "FETCH_HEAD"], cwd=ROOT, check=True,
        stdout=subprocess.PIPE,
    ).stdout
    with zipfile.ZipFile(io.BytesIO(archive)) as zf:
        zf.extractall(DEMOS)
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
    env = dict(os.environ, TIGRE_DEMOS_RUN="1")
    cmd = [sys.executable, "-m", "sphinx", "-b", "html", str(DOCS), str(DOCS / "_build" / "html")]
    print(" ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)


THUMB_WIDTH = 230  # thumbnails are shown at about this width on the section pages


def shrink_thumbnails() -> None:
    """Make the animated thumbnails light enough for the section pages.

    sphinx-gallery uses a demo's animation, at full size, as its thumbnail
    (up to 1 MB each). Keep the animation but make it smaller, with at most 32
    frames and 64 colours, and start it half-way through: the first slice of a
    volume is often empty, which made some thumbnails look black. Thumbnails
    reused from the previous render are already small and are left alone.
    """
    from PIL import Image, ImageSequence

    for path in DEMOS.glob("*/images/thumb/*.gif"):
        im = Image.open(path)
        if im.width <= THUMB_WIDTH:
            continue
        frames = [f.convert("RGB") for f in ImageSequence.Iterator(im)]
        duration = im.info.get("duration", 100)
        if len(frames) > 1:
            frames = frames[len(frames) // 2:] + frames[: len(frames) // 2]
            step = -(-len(frames) // 32)
            frames = frames[::step]
            duration *= step
        height = round(frames[0].height * THUMB_WIDTH / frames[0].width)
        frames = [f.resize((THUMB_WIDTH, height), Image.LANCZOS) for f in frames]
        palette = frames[0].quantize(colors=64, method=Image.Quantize.MEDIANCUT)
        frames = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
        frames[0].save(path, save_all=True, append_images=frames[1:], loop=0,
                       duration=duration, optimize=True, disposal=1)


def commit() -> str:
    """Commit docs/demos as the new tip of the docs-demos branch."""
    (DEMOS / ".gitattributes").write_text("* -text\n", encoding="utf-8")
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ, GIT_INDEX_FILE=str(pathlib.Path(tmp) / "index"))
        git("--work-tree", str(DEMOS), "add", "--all", "--force", ".", env=env)
        tree = git("write-tree", env=env)
    # One commit with no history: each render replaces the last, so old images
    # do not pile up in the repository. Publishing therefore needs --force.
    head = git("rev-parse", "--short", "HEAD")
    sha = git("commit-tree", tree, "-m", f"Rendered demos from {head}")
    git("update-ref", f"refs/heads/{BRANCH}", sha)
    return sha


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--fresh", action="store_true", help="re-run every demo")
    parser.add_argument("--push", action="store_true", help="push and redeploy (needs gh)")
    args = parser.parse_args()

    check_tigre()
    if DEMOS.exists():
        shutil.rmtree(DEMOS)
    if not args.fresh:
        restore_previous()
    build()
    shrink_thumbnails()
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
