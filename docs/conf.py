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
    "sphinxext.opengraph",
    "sphinx_sitemap",
    "notfound.extension",
]

# -- API reference -------------------------------------------------------------
# autodoc imports tigre without building it: the compiled CUDA extensions and
# optional dependencies are replaced by mocks.
autodoc_mock_imports = [
    *("tigre.utilities.cuda_interface." + m for m in (
        "_Ax", "_Atb", "_AwminTV", "_minTV", "_minPICCS", "_tv_proximal", "_gpuUtils", "_RandomNumberGenerator")),
    "h5py", "tqdm", "torch",
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

# Wheels attached to each GitHub release: their platform tag comes from the
# cibuildwheel image (manylinux_2_28 since #827). If a release changes it, update
# _WHEEL_PLATFORM, or the install commands on the Python install page will 404.
_WHEEL_PYTHONS = ["311", "312", "313", "314"]
_WHEEL_PLATFORM = "manylinux_2_24_x86_64.manylinux_2_28_x86_64"


def _wheel_commands():
    """Write the copy-paste pip commands for the latest release's wheels."""
    base = f"https://github.com/CERN/TIGRE/releases/download/v{latest_version}"
    lines = ["::::{tab-set}", ":sync-group: python", ""]
    for py in _WHEEL_PYTHONS:
        tag = f"cp{py}"
        wheel = f"tigre-{latest_version}-{tag}-{tag}-{_WHEEL_PLATFORM}.whl"
        lines += [f":::{{tab-item}} Python 3.{py[1:]}", f":sync: py3{py[1:]}", "",
                  "```sh", f"pip install {base}/{wheel}", "```", ":::", ""]
    lines.append("::::")
    out = pathlib.Path(__file__).resolve().parent / "_generated" / "wheels.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


_wheel_commands()
_n_papers = len(re.findall(r"^- ", (pathlib.Path(__file__).resolve().parent / "papers.yaml")
                            .read_text(encoding="utf-8"), flags=re.M))
myst_substitutions = {
    "release": release,
    "latest_version": latest_version,
    "default_mode": "auto",  # light or dark, following the reader's system setting
    "latest_date": latest_date,
    "papers_count": f"{_n_papers // 10 * 10}+",
    "latest_download": f"[download it](https://github.com/CERN/TIGRE/releases/tag/v{latest_version})",
}

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "_generated", "demos/**/*.ipynb", "demos/**/*.py.md5", "papers.yaml", "Thumbs.db", ".DS_Store", "README.md",
                    "technical-notes/**/README.md"]

nitpicky = True
linkcheck_anchors = False
linkcheck_ignore = [
    r"https://doi\.org/.*",  # publishers often block automated link checkers
    r"https://iopscience\.iop\.org/.*",
    r"https://eprints\.soton\.ac\.uk/.*",  # refuses automated requests (401/403); pages are public
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
    # Thumbnail for demos shown without output (they need data or hardware the
    # website build does not have): a grey TIGRE icon.
    "default_thumb_file": str(pathlib.Path(__file__).parent / "_static" / "demo-no-output.png"),
    # Notebooks open in Google Colab (see _colab_button below): their first cell
    # installs the latest release's wheel there, and does nothing elsewhere.
    "first_notebook_cell": (
        "# On Google Colab (with a GPU runtime), install TIGRE's prebuilt wheel.\n"
        "# Elsewhere this cell does nothing: install TIGRE as described at\n"
        "# https://cern.github.io/TIGRE/install/index.html\n"
        "import sys, subprocess\n"
        "if 'google.colab' in sys.modules:\n"
        "    _py = f'cp{sys.version_info.major}{sys.version_info.minor}'\n"
        "    _wheel = (f'https://github.com/CERN/TIGRE/releases/download/v" + latest_version + "/'\n"
        "              f'tigre-" + latest_version + "-{_py}-{_py}-" + _WHEEL_PLATFORM + ".whl')\n"
        "    subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', _wheel], check=True)"
    ),
    # do not print the value of a demo's last expression (e.g. "<plotImg object at 0x...>")
    "capture_repr": (),
}

# -- HTML output -------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = "TIGRE"
html_static_path = ["_static"]
html_css_files = ["custom.css", "landing.css"]
html_js_files = ["compare.js"]
html_baseurl = "https://cern.github.io/TIGRE/"
html_favicon = "_static/logo/tigre-icon.svg"
html_extra_path = ["_extra"]  # robots.txt

# Link previews (Slack, Teams, social media, email): title, description and image.
ogp_site_url = "https://cern.github.io/TIGRE/"
ogp_site_name = "TIGRE"
ogp_image = "https://cern.github.io/TIGRE/_static/logo/tigre-social-preview.png"
ogp_image_alt = "TIGRE: fast, flexible 3D CT reconstruction on GPUs. A coffee bean reconstructed with FDK and CGLS, split by an orange divider."
ogp_social_cards = {"enable": False}

# sitemap.xml for search engines (from html_baseurl).
sitemap_url_scheme = "{link}"
sitemap_excludes = ["search.html", "genindex.html", "404.html"]

# A 404 page in the site's style. GitHub Pages serves 404.html for any missing
# page, from any folder, so its links must be absolute.
notfound_urls_prefix = "/TIGRE/"
notfound_context = {
    "title": "Page not found",
    "body": (
        "<h1>Page not found</h1>"
        "<p>This page does not exist, or it has moved. The website was rebuilt in "
        "November 2026, so old links may point to pages that are now elsewhere.</p>"
        "<p>Try the search box above, or start from the <a href='/TIGRE/'>home page</a>, "
        "the <a href='/TIGRE/install/index.html'>install guide</a> or the "
        "<a href='/TIGRE/user-guide/index.html'>user guide</a>.</p>"
    ),
}
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
    "footer_center": ["footer-links"],
    "footer_end": [],
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
    "default_mode": "auto",  # light or dark, following the reader's system setting
    "github_user": "CERN",
    "github_repo": "TIGRE",
    "github_version": "master",
    "doc_path": "docs",
}

# Single-page sections have nothing to show in the left sidebar.
html_sidebars = {
    page: []
    for page in [
        "index", "gallery/index", "citing", "papers", "ecosystem", "roadmap", "help", "contributing",
        "changelog", "acknowledgements",
    ]
}


# -- Open in Colab --------------------------------------------------------------
# Each Python demo page gets a button that opens its notebook in Google Colab.
# The notebooks are on the docs-demos branch next to the rendered pages.
_COLAB_SKIP = {"d21_ScannerDataLoader", "d24_CurvedDetector"}  # need your own data


def _colab_button(app, docname, source):
    parts = docname.split("/")
    if len(parts) != 3 or parts[0] != "demos" or parts[2] in _COLAB_SKIP:
        return
    notebook = pathlib.Path(app.srcdir) / f"{docname}.ipynb"
    if not notebook.exists():
        return  # MATLAB-only demo
    url = (f"https://colab.research.google.com/github/CERN/TIGRE/blob/docs-demos/"
           f"{parts[1]}/{parts[2]}.ipynb")
    button = (
        f"\n.. button-link:: {url}\n"
        "   :color: secondary\n"
        "   :outline:\n\n"
        "   :octicon:`play` Open in Colab\n\n"
        "Runs in your browser on Google Colab: choose a GPU runtime "
        "(*Runtime > Change runtime type*), then run the cells in order.\n"
    )
    text = source[0]
    marker = text.find("\nSource: ")
    if marker == -1:
        return
    end = text.find("\n", marker + 1)
    source[0] = text[:end + 1] + button + text[end + 1:]


# -- Demo thumbnails -------------------------------------------------------------
# sphinx-gallery writes the thumbnails with an empty alt text. Give each one the
# title of its demo, so screen readers and image search know what it shows.
_THUMB = re.compile(
    r'(<div class="sphx-glr-thumbcontainer"[^>]*>)<img alt=""'
    r'(.*?<div class="sphx-glr-thumbnail-title">)(.*?)(</div>)',
    re.S,
)


def _thumbnail_alt(app, exception):
    if exception is not None or app.builder.format != "html":
        return
    for page in pathlib.Path(app.outdir, "demos").rglob("index.html"):
        text = page.read_text(encoding="utf-8")
        new = _THUMB.sub(
            lambda m: f'{m[1]}<img alt="{m[3].strip().replace(chr(34), "&quot;")}"{m[2]}{m[3]}{m[4]}', text
        )
        if new != text:
            page.write_text(new, encoding="utf-8")


def _demo_thumbnails(app, exception):
    """Tidy the demo thumbnails on the demos index.

    - Animated thumbnails become still. Demos that animate a volume get an
      animated GIF thumbnail, so the index had a dozen images moving at once;
      the animations stay on the demo pages. The kept frame is the one with
      most mid-grey pixels, which shows most of the object (the first and last
      slices are often empty).
    - Demos shown without output get a grey TIGRE icon instead of
      sphinx-gallery's blue and yellow placeholder. (``default_thumb_file``
      does the same for newly rendered demos; this covers the stored ones.)
    """
    if exception is not None or app.builder.format != "html":
        return
    import numpy as np
    from PIL import Image

    images = pathlib.Path(app.outdir, "_images")
    placeholder = pathlib.Path(app.confdir, "_static", "demo-no-output.png")

    def mid_grey(frame):
        grey = np.asarray(frame.convert("L"))
        return int(((grey > 40) & (grey < 220)).sum())

    for thumb in images.glob("sphx_glr_*_thumb.gif"):
        with Image.open(thumb) as im:
            frames = getattr(im, "n_frames", 1)
            if frames < 2:
                continue
            scores = []
            for k in range(frames):
                im.seek(k)
                scores.append(mid_grey(im))
            im.seek(int(np.argmax(scores)))
            still = im.convert("RGB")
        still.save(thumb)

    for thumb in images.glob("sphx_glr_*_thumb.png"):
        with Image.open(thumb) as im:
            rgb = np.asarray(im.convert("RGB")).astype(int)
            size = im.size
        blue = (np.abs(rgb - [173, 216, 230]).sum(-1) < 40).mean()
        yellow = (np.abs(rgb - [255, 228, 181]).sum(-1) < 40).mean()
        if blue > 0.005 and yellow > 0.1:
            with Image.open(placeholder) as ph:
                ph.convert("RGB").resize(size).save(thumb)


def setup(app):
    app.connect("source-read", _colab_button)
    app.connect("build-finished", _thumbnail_alt)
    app.connect("build-finished", _demo_thumbnails)
