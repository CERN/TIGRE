"""Sphinx configuration for the TIGRE website and documentation."""

import datetime
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
exclude_patterns = ["_build", "_generated", "papers.yaml", "Thumbs.db", ".DS_Store", "README.md"]

nitpicky = True
linkcheck_anchors = False
linkcheck_ignore = [
    r"https://doi\.org/.*",  # publishers often block automated link checkers
    r"https://iopscience\.iop\.org/.*",
]

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
