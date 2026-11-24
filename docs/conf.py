"""Sphinx configuration for the TIGRE website and documentation."""

import datetime
import os
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "_ext"))
# the tigre package, for the API reference (its CUDA modules are mocked below)
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "Python"))

# -- Project information -----------------------------------------------------

project = "TIGRE"
author = "Ander Biguri and the TIGRE contributors"
copyright = f"2015-{datetime.date.today().year}, University of Bath and CERN (BSD licence)"


def _version_from_cff():
    """Read the release version from CITATION.cff (single source of truth)."""
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r'^version:\s*"?([^"\n]+)"?', text, flags=re.M)
    return match.group(1).strip() if match else "dev"


release = _version_from_cff()
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
    "tigre_citations",
    "tigre_papers",
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "tigre_matlab_api",
]

# -- API reference -------------------------------------------------------------
# autodoc imports tigre without building it: the compiled CUDA extensions and
# optional dependencies are replaced by mocks.
autodoc_mock_imports = [
    "_Ax", "_Atb", "_AwminTV", "_minTV", "_minPICCS", "_tv_proximal", "_tvdenoising",
    "_gpuUtils", "_RandomNumberGenerator", "h5py", "tqdm", "torch",
]
autodoc_typehints = "none"
autodoc_member_order = "bysource"
autosummary_generate = False
napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_use_rtype = False
# Docstring types ("numpy.ndarray", "optional", {"a", "b"}) are descriptions,
# not cross-references; keep nitpicky mode for everything else.
nitpick_ignore_regex = [("py:class", r".*")]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "substitution",
]
myst_heading_anchors = 3


def _changelog():
    """Copy changelog.md for the site, with #123 turned into links to GitHub.

    Returns the newest released version and its date, for the landing page.
    """
    text = (ROOT / "changelog.md").read_text(encoding="utf-8")
    body = text[text.index("\n## "):]  # the intro is for contributors
    body = re.sub(r"## Unreleased\n+(?=## )", "", body)  # nothing unreleased yet
    linked = re.sub(r"(?<![\w/\[`])#(\d{2,4})\b",
                    r"[#\1](https://github.com/CERN/TIGRE/pull/\1)", body)
    out = pathlib.Path(__file__).resolve().parent / "_generated" / "changelog.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(linked, encoding="utf-8")
    match = re.search(r"^## (\d+\.\d+\.\d+) \((\d{4}-\d{2}-\d{2})\)", text, flags=re.M)
    if not match:
        return release, ""
    day = datetime.date.fromisoformat(match.group(2))
    return match.group(1), f"{day.day} {day:%B %Y}"


latest_version, latest_date = _changelog()
_n_papers = len(re.findall(r"^- ", (pathlib.Path(__file__).resolve().parent / "papers.yaml")
                            .read_text(encoding="utf-8"), flags=re.M))
myst_substitutions = {
    "release": release,
    "latest_version": latest_version,
    "latest_date": latest_date,
    "papers_count": f"{_n_papers // 10 * 10}+",
    "latest_download": f"[Download](https://github.com/CERN/TIGRE/releases/tag/v{latest_version})",
}

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "_generated", "demos/**/*.ipynb", "demos/**/*.py.md5", "papers.yaml", "Thumbs.db", ".DS_Store", "README.md",
                    "technical-notes/**/README.md"]

nitpicky = True
linkcheck_anchors = False
linkcheck_ignore = [
    r"https://doi\.org/.*",  # publishers often block automated link checkers
    r"https://iopscience\.iop\.org/.*",
]

# -- Demos (sphinx-gallery) ------------------------------------------------
# The demo pages are generated from Python/demos and MATLAB/Demos (see
# _ext/tigre_demos.py). Running the demos needs TIGRE and a GPU, so a normal
# build only shows the code. Ander renders the outputs on a GPU machine with
# docs/tools/build_demos.py and keeps them on the docs-demos branch; the
# deploy copies that branch into docs/demos, and sphinx-gallery reuses each
# rendered page as long as its demo has not changed.

import tigre_demos  # noqa: E402

_demos_src = tigre_demos.generate(pathlib.Path(__file__).resolve().parent / "_generated" / "demos_src")

sphinx_gallery_conf = {
    "examples_dirs": [str(_demos_src)],
    "gallery_dirs": ["demos"],
    "filename_pattern": tigre_demos.run_pattern(),
    "plot_gallery": os.environ.get("TIGRE_DEMOS_RUN", "") == "1",
    "example_extensions": {".py", ".m"},
    "filetype_parsers": {".m": "Matlab"},
    "within_subsection_order": "tigre_demos.CatalogueOrder",
    "download_all_examples": False,
    "remove_config_comments": True,
    "show_memory": False,
    "write_computation_times": False,
    "image_scrapers": ("tigre_demos.scraper",),
    "first_notebook_cell": None,
    # do not print the value of a demo's last expression (e.g. "<plotImg object at 0x...>")
    "capture_repr": (),
}

# -- HTML output -------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = "TIGRE"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_baseurl = "https://cern.github.io/TIGRE/"
html_favicon = "_static/logo/tigre-icon.svg"
templates_path = ["_templates"]

html_theme_options = {
    "github_url": "https://github.com/CERN/TIGRE",
    "logo": {
        "image_light": "_static/logo/tigre-logo.svg",
        "image_dark": "_static/logo/tigre-logo-dark.svg",
        "alt_text": "TIGRE - home",
    },
    "header_links_before_dropdown": 6,
    "navbar_align": "left",
    "navbar_end": ["navbar-version", "theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "use_edit_page_button": True,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "icon_links": [
        {
            "name": "Discussions",
            "url": "https://github.com/CERN/TIGRE/discussions",
            "icon": "fa-solid fa-comments",
        },
    ],
}

html_context = {
    "latest_version": latest_version,
    "github_user": "CERN",
    "github_repo": "TIGRE",
    "github_version": "master",
    "doc_path": "docs",
}

# Single-page sections have nothing to show in the left sidebar.
html_sidebars = {
    page: [] for page in ["index", "citing", "papers", "changelog", "acknowledgements"]
}
