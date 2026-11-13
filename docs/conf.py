"""Sphinx configuration for the TIGRE website and documentation."""

import datetime
import os
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "_ext"))

# -- Project information -----------------------------------------------------

project = "TIGRE"
author = "Ander Biguri and the TIGRE contributors"
copyright = f"2015-{datetime.date.today().year}, University of Bath and CERN"


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
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "substitution",
]
myst_heading_anchors = 3
myst_substitutions = {"release": release}

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "_generated", "gallery/**/*.ipynb", "gallery/**/*.py.md5", "papers.yaml", "Thumbs.db", ".DS_Store", "README.md"]

nitpicky = True
linkcheck_anchors = False
linkcheck_ignore = [
    r"https://doi\.org/.*",  # publishers often block automated link checkers
    r"https://iopscience\.iop\.org/.*",
]

# -- Gallery -----------------------------------------------------------------
# The gallery is generated from Python/demos and MATLAB/Demos (see
# _ext/tigre_gallery.py). Running the demos needs TIGRE and a GPU, so a normal
# build only shows the code. Ander renders the outputs on a GPU machine with
# docs/tools/build_gallery.py and keeps them on the docs-gallery branch; the
# deploy copies that branch into docs/gallery, and sphinx-gallery reuses each
# rendered page as long as its demo has not changed.

import tigre_gallery  # noqa: E402

_gallery_src = tigre_gallery.generate(pathlib.Path(__file__).resolve().parent / "_generated" / "gallery_src")

sphinx_gallery_conf = {
    "examples_dirs": [str(_gallery_src)],
    "gallery_dirs": ["gallery"],
    "filename_pattern": tigre_gallery.run_pattern(),
    "plot_gallery": os.environ.get("TIGRE_GALLERY_RUN", "") == "1",
    "example_extensions": {".py", ".m"},
    "filetype_parsers": {".m": "Matlab"},
    "matplotlib_animations": (True, "jshtml"),
    "within_subsection_order": "tigre_gallery.CatalogueOrder",
    "download_all_examples": False,
    "remove_config_comments": True,
    "show_memory": False,
    "write_computation_times": False,
    "image_scrapers": ("tigre_gallery.scraper",),
    "first_notebook_cell": None,
}

# -- HTML output -------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = "TIGRE"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_baseurl = "https://cern.github.io/TIGRE/"

html_theme_options = {
    "github_url": "https://github.com/CERN/TIGRE",
    "header_links_before_dropdown": 6,
    "navbar_align": "left",
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
    "github_user": "CERN",
    "github_repo": "TIGRE",
    "github_version": "master",
    "doc_path": "docs",
}

# Single-page sections have nothing to show in the left sidebar.
html_sidebars = {
    page: [] for page in ["index", "citing", "papers"]
}
