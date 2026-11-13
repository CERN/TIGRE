"""Turn TIGRE's demos into sources for sphinx-gallery.

The demos in ``Python/demos`` and ``MATLAB/Demos`` are meant to be run as they
are, so they are not edited for the website. Instead, ``generate()`` writes a
converted copy of each one to ``docs/_generated/gallery_src``:

* the demo's own header (title, license block) is replaced by the title and
  description in ``CATALOGUE`` below, as sphinx-gallery expects;
* ``#%% Step`` cell markers become sphinx-gallery text cells, so each step gets
  a heading on the page; other comments stay in the code;
* a Python demo with a MATLAB twin gets a drop-down with the MATLAB file;
* MATLAB-only demos are shown as code (they are not run).

``conf.py`` calls ``generate()`` before sphinx-gallery reads the sources. The
output only depends on the demo files and this catalogue, so a rendered
gallery (``docs/gallery``, kept on the ``docs-gallery`` branch) stays valid
until one of them changes; see ``docs/tools/build_gallery.py``.
"""
from __future__ import annotations

import pathlib
import re
import shutil
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[2]
PY_DEMOS = ROOT / "Python" / "demos"
M_DEMOS = ROOT / "MATLAB" / "Demos"

# Demos that need data, hardware or packages a gallery build does not have
# (d24 is a template: it needs your own curved-detector projections).
NOT_RUN = {"d21_ScannerDataLoader.py", "d24_CurvedDetector.py", "d25_Pytorch.py"}

# (folder, heading, intro, [(python file or None, MATLAB file or None, title, description)])
CATALOGUE = [
    (
        "01_basics",
        "Getting started",
        "Describe a scanner, simulate data, reconstruct and look at the result.",
        [
            ("d01_CreateGeometry.py", "d01_CreateGeometry.m", "Describe the scanner geometry",
             "Create a cone-beam geometry field by field, and the default geometries "
             "that save typing. See :doc:`/user-guide/geometry`."),
            ("d02_SampleData.py", "d02_SampleData.m", "Load the sample data",
             "Load the head phantom and the other images that ship with TIGRE."),
            ("d03_generateData.py", "d03_generateData.m", "Simulate projections and noise",
             "Forward-project an image and add realistic Poisson and Gaussian noise."),
            ("d04_SimpleReconstruction.py", "d04_SimpleReconstruction.m", "A first reconstruction",
             "Reconstruct simulated data with FDK and OS-SART and compare the results."),
            ("d05_plottingFunctions.py", "d05_plottingFunctions.m", "Plot images and projections",
             "All the options of ``plotimg`` and ``plotproj``."),
            ("d10_MeasureQuality.py", "d10_QualityMeasurements.m", "Measure image quality",
             "Compute error and similarity measures against a reference image, also "
             "during the iterations of an algorithm."),
            ("d11_PostProcessing.py", "d11_PostProcessing.m", "Post-process a reconstruction",
             "Denoise a reconstruction with total variation after it is computed."),
        ],
    ),
    (
        "02_algorithms",
        "Algorithms",
        "Each family of reconstruction algorithms on the same data. "
        "See :doc:`/user-guide/algorithms` for how to choose.",
        [
            ("d06_Algorithms01.py", "d06_Algorithms01.m", "FDK and its filters",
             "Filtered backprojection for cone beam, with the available filters."),
            ("d07_Algorithms02.py", "d07_Algorithms02.m", "SART, OS-SART and SIRT",
             "The SART family, the options they share, and how they compare."),
            ("d08_Algorithms03.py", "d08_Algorithms03.m", "Krylov subspace methods",
             "CGLS, LSQR, LSMR, hybrid LSQR and the GMRES variants."),
            ("d09_Algorithms04.py", "d09_Algorithms04.m", "Total variation algorithms",
             "ASD-POCS and its variants, and their parameters."),
            ("d20_Algorithms05.py", "d20_Algorithms05.m", "FISTA with a TV proximal step",
             "FISTA, which reaches a TV-regularised solution faster than gradient methods."),
            (None, "d24_Algorithms06.m", "Prior images with PICCS",
             "Use a prior image of the same object, such as a planning CT."),
        ],
    ),
    (
        "03_geometry",
        "Geometries and operators",
        "Projection operators and the scanner geometries TIGRE can describe.",
        [
            ("d12_ProjectionOperations.py", "d12_ProjectionOperators.m", "Projectors and backprojectors",
             "Siddon and interpolated projection, FDK-weighted and matched backprojection. "
             "See :doc:`/user-guide/projectors`."),
            ("d13_HelicalGeometry.py", "d13_HelicalGeometry.m", "Helical scans",
             "Move the image along the rotation axis at each angle."),
            ("d14_Offsets.py", "d14_Offsets.m", "Image and detector offsets",
             "Shift the image and the detector, per scan or per angle."),
            ("d15_3Dparallel.py", "d15_3Dparallel.m", "Parallel-beam geometry",
             "Reconstruct 3D parallel-beam data, as from a synchrotron."),
            ("d16_2Dtomography.py", "d16_2Dtomography.m", "2D tomography",
             "Fan-beam and parallel-beam 2D problems, as one-slice 3D problems."),
            ("d17_DetectorRotation.py", "d17_DetectorRotation.m", "Detector rotation",
             "Describe a detector that is rotated or tilted."),
            ("d18_ArbitraryAxisOfRotation.py", "d18_ArbitraryAxisOfRotation.m",
             "Arbitrary axis of rotation",
             "Trajectories that are not a circle around z, given as Euler angles."),
            ("d19_DBTexample.py", "d19_DBTexample.m", "Digital breast tomosynthesis",
             "A limited-angle cone-beam scan, reconstructed with FDK and OS-SART."),
            ("d24_CurvedDetector.py", "d024_CurvedDetector.m", "Curved detectors",
             "Resample curved-detector data onto a flat detector, then reconstruct."),
            ("d26_ShortScanWangWeights.py", "d26_ShortScanWangWeights.m",
             "Shifted detector on a short scan",
             "Why Wang weights are skipped on short scans, shown on a uniform cylinder."),
        ],
    ),
    (
        "04_data",
        "Real data, GPUs and other tools",
        "Loading scanner data, choosing GPUs and using TIGRE from other software. "
        "These demos need data or hardware that the website build does not have, "
        "so they are shown without output.",
        [
            ("d21_ScannerDataLoader.py", "d21_ScannerDataLoader.m", "Load scanner data",
             "Read Nikon, Bruker, YXLON, Diondo and Varian datasets. "
             "See :doc:`/user-guide/loading-data`."),
            ("d22_ListGPUs.py", "d22_ListGPUs.m", "List and choose GPUs",
             "Find the installed GPUs and choose which ones TIGRE uses. "
             "See :doc:`/user-guide/gpus`."),
            ("d25_Pytorch.py", None, "Use TIGRE from PyTorch",
             "Differentiable TIGRE operators inside a learned reconstruction network."),
            (None, "d23_DualEnergy.m", "Dual-energy CBCT",
             "Material decomposition from two scans at different tube voltages."),
            (None, "SophiaBeads_demo.m", "Real micro-CT data: SophiaBeads",
             "Reconstruct the SophiaBeads dataset, a tube of glass beads."),
            (None, "pRad_demo.m", "Proton radiography",
             "Compute an optimised proton radiography. See :doc:`/user-guide/proton-ct`."),
        ],
    ),
]

_PY_CELL = re.compile(r"^(?:#+\s*%%|##%)\s*(.*)$")
_M_CELL = re.compile(r"^%%\s*(.*)$")


def _rst_escape(text: str) -> str:
    return re.sub(r"([\\`*_|])", r"\\\1", text.strip())


def _heading(title: str, char: str = "=") -> str:
    return f"{title}\n{char * len(title)}"


def _read(path: pathlib.Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").split("\n")


def _strip_header(lines: list[str], comment: str) -> list[str]:
    """Drop a leading docstring and the leading block of comments and blank lines."""
    i = 0
    if comment == "#" and lines and lines[0].lstrip().startswith(('"""', "'''")):
        quote = lines[0].lstrip()[:3]
        if lines[0].count(quote) >= 2 and len(lines[0].strip()) > 3:
            i = 1
        else:
            i = 1
            while i < len(lines) and quote not in lines[i]:
                i += 1
            i += 1
    cell = _PY_CELL if comment == "#" else _M_CELL
    start = i
    while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith(comment)):
        # a step marker after the demo's title line starts the content
        if i > start and cell.match(lines[i]):
            break
        i += 1
    return lines[i:]


def _convert_body(lines: list[str], language: str) -> str:
    cell = _PY_CELL if language == "python" else _M_CELL
    marker = "# %%" if language == "python" else "%%"
    prefix = "# " if language == "python" else "% "
    out = []
    for line in lines:
        match = cell.match(line)
        if match:
            title = match.group(1).strip().rstrip(":").strip()
            out.append("")
            if title:
                if language == "python":
                    # a sub-heading, like MATLAB's "%%" sections
                    heading = _rst_escape(title)
                    out += [marker, prefix + heading, prefix + "-" * len(heading)]
                else:
                    out.append(f"{marker} {_rst_escape(title)}")
            else:
                out.append(marker)
            # a blank line ends the text cell; the demo's own comments stay as code
            out.append("")
        else:
            out.append(line.rstrip())
    text = "\n".join(out).strip("\n")
    return re.sub(r"\n{3,}", "\n\n", text) + "\n"


def _matlab_dropdown(m_file: str) -> str:
    return textwrap.dedent(
        f"""
        .. dropdown:: MATLAB version of this demo
           :icon: code

           .. literalinclude:: /../MATLAB/Demos/{m_file}
              :language: matlab
        """
    )


def _python_example(py: str, m: str | None, title: str, description: str) -> str:
    doc = f"{_heading(title)}\n\n{description}\n"
    doc += f"\nSource: `Python/demos/{py} <https://github.com/CERN/TIGRE/blob/master/Python/demos/{py}>`__"
    if m:
        doc += f" · `MATLAB/Demos/{m} <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/{m}>`__\n"
        doc += _matlab_dropdown(m)
    else:
        doc += "\n\nThis demo is only available in Python.\n"
    if py in NOT_RUN:
        doc += "\nThis demo is not run when the website is built.\n"
    body = _convert_body(_strip_header(_read(PY_DEMOS / py), "#"), "python")
    return f'"""\n{doc.rstrip()}\n"""\n\n{body}'


def _matlab_example(m: str, title: str, description: str) -> str:
    doc = (
        f"{_heading(title)}\n\n{description}\n\n"
        f"Source: `MATLAB/Demos/{m} <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/{m}>`__\n\n"
        "This demo is only available in MATLAB, so it is shown without output."
    )
    header = "\n".join("% " + line if line else "%" for line in doc.split("\n"))
    body = _convert_body(_strip_header(_read(M_DEMOS / m), "%"), "matlab")
    return f"{header}\n\n{body}"


def generate(out_dir: pathlib.Path) -> pathlib.Path:
    """Write the gallery sources to out_dir and return it."""
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    (out_dir / "GALLERY_HEADER.rst").write_text(
        _heading("Gallery")
        + "\n\nWorked examples, generated from TIGRE's demos. Each page shows the Python "
        "demo and its output; open the drop-down at the top of a page for the MATLAB "
        "version. The demos themselves are in "
        "`Python/demos <https://github.com/CERN/TIGRE/tree/master/Python/demos>`__ and "
        "`MATLAB/Demos <https://github.com/CERN/TIGRE/tree/master/MATLAB/Demos>`__.\n",
        encoding="utf-8",
        newline="\n",
    )
    for folder, heading, intro, demos in CATALOGUE:
        section = out_dir / folder
        section.mkdir()
        (section / "GALLERY_HEADER.rst").write_text(
            f"{_heading(heading, '-')}\n\n{intro}\n", encoding="utf-8", newline="\n"
        )
        for py, m, title, description in demos:
            if py:
                (section / py).write_text(
                    _python_example(py, m, title, description), encoding="utf-8", newline="\n"
                )
            else:
                (section / m).write_text(
                    _matlab_example(m, title, description), encoding="utf-8", newline="\n"
                )
    return out_dir


def run_pattern() -> str:
    """Regular expression for the examples sphinx-gallery should execute."""
    names = [py for _, _, _, demos in CATALOGUE for py, *_ in demos if py and py not in NOT_RUN]
    return "(" + "|".join(re.escape(n) for n in names) + ")$"


class CatalogueOrder:
    """Sort key for sphinx-gallery: the order of CATALOGUE, not file names."""

    def __init__(self, src_dir):
        self.src_dir = src_dir
        self.order = {}
        for _, _, _, demos in CATALOGUE:
            for i, (py, m, *_) in enumerate(demos):
                self.order[py or m] = i

    def __call__(self, filename):
        return self.order.get(pathlib.Path(filename).name, len(self.order))

    def __repr__(self):
        return f"<{self.__class__.__name__}>"


def scraper(block, block_vars, gallery_conf, **kwargs):
    """sphinx-gallery image scraper that shows TIGRE's plot animations as GIFs.

    ``tigre.plotimg`` and ``tigre.plotproj`` animate through the slices of a
    volume, and attach the animation to its figure. sphinx-gallery's own
    animation support embeds every frame in the page as HTML, which made some
    pages over 30 MB, so save those animations as GIF images instead and let
    the standard matplotlib scraper handle the remaining (still) figures.
    """
    from pathlib import PurePosixPath

    import matplotlib.pyplot as plt
    from sphinx_gallery.scrapers import figure_rst, matplotlib_scraper

    paths = block_vars["image_path_iterator"]
    gifs = []
    for num in plt.get_fignums():
        fig = plt.figure(num)
        anim = getattr(fig, "_tigre_animation", None)
        if anim is None:
            continue
        gif = PurePosixPath(next(paths)).with_suffix(".gif")
        interval = getattr(anim, "_interval", 200) or 200
        anim.save(str(gif), writer="pillow", fps=max(1, round(1000 / interval)), dpi=72)
        gifs.append(gif)
        plt.close(fig)
    rst = figure_rst(gifs, gallery_conf["src_dir"]) if gifs else ""
    return rst + matplotlib_scraper(block, block_vars, gallery_conf, **kwargs)
