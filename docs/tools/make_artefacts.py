"""Make the images of the website's "Why does my reconstruction look wrong?" page.

Each case simulates projections of the head phantom with one common mistake (a
wrong centre of rotation, a missing log, a truncated detector, ...) and saves a
figure with the correct reconstruction next to the wrong one, as
``<case>.png`` in the output folder. It needs TIGRE and a GPU; it takes a few
minutes.

    python docs/tools/make_artefacts.py --out <folder>

Run it from outside TIGRE's source folder, or with the repository root as the
current folder (``python docs/tools/make_artefacts.py``): it imports the installed
TIGRE either way. The images go to ``docs/user-guide/images/artefacts``.
"""

import argparse
import copy
import pathlib
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import tigre  # noqa: E402
import tigre.algorithms as algs  # noqa: E402
from tigre.utilities import CTnoise, sample_loader  # noqa: E402

N = 256  # image size (voxels per side)
WINDOW = (0.0, 1.0)  # display window for the head phantom


def make_geometry():
    """A cone-beam scanner whose detector sees the whole 256 mm image."""
    geo = tigre.geometry_default(high_resolution=False)
    geo.DSD = 1536.0
    geo.DSO = 1000.0
    geo.nVoxel = np.array([N, N, N])
    geo.sVoxel = np.array([256.0, 256.0, 256.0])
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    geo.nDetector = np.array([N, N])
    geo.dDetector = np.array([1.6, 1.6])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offOrigin = np.array([0.0, 0.0, 0.0])
    geo.offDetector = np.array([0.0, 0.0])
    geo.rotDetector = np.array([0.0, 0.0, 0.0])
    geo.COR = 0.0
    return geo


def axial(img):
    return img[img.shape[0] // 2]


def coronal(img):
    return img[:, img.shape[1] // 2, :]


def figure(out, name, panels, window=WINDOW):
    """Save panels [(title, 2D image, window or None)] side by side."""
    fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 3.5))
    for ax, (title, img, win) in zip(np.atleast_1d(axes), panels):
        lo, hi = win if win is not None else window
        ax.imshow(img, cmap="gray", vmin=lo, vmax=hi)
        ax.set_title(title, fontsize=11)
        ax.axis("off")
    fig.tight_layout()
    path = out / f"{name}.png"
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def diff_window(img):
    m = float(np.percentile(np.abs(img), 99.5))
    return -m, m


def auto_window(img):
    lo, hi = np.percentile(img, [1, 99])
    return float(lo), float(hi)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="docs/user-guide/images/artefacts", type=pathlib.Path)
    parser.add_argument("--only", nargs="*", help="make only these cases")
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    only = set(args.only or [])

    def wanted(name):
        return not only or name in only

    t0 = time.time()
    base = make_geometry()

    def G():
        # In TIGRE 3.3.1 and earlier, the iterative algorithms wrote per-angle copies
        # of the geometry fields into the geo they were given, so every call gets a copy.
        return copy.deepcopy(base)

    geo = base
    head = sample_loader.load_head_phantom(geo.nVoxel).astype(np.float32)
    angles = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    proj = tigre.Ax(head, G(), angles)
    good = algs.fdk(proj, G(), angles)
    print(f"reference done ({time.time() - t0:.0f} s)")

    # 1. Centre of rotation: the data have a COR shift that the reconstruction ignores.
    if wanted("cor"):
        g = G()
        g.COR = 1.5  # mm, about one voxel
        p = tigre.Ax(head, g, angles)
        g.COR = 4.0
        p4 = tigre.Ax(head, g, angles)
        figure(out, "cor", [
            ("Correct", axial(good), None),
            ("COR off by 1.5 mm", axial(algs.fdk(p, G(), angles)), None),
            ("COR off by 4 mm", axial(algs.fdk(p4, G(), angles)), None),
        ])

    # 2. Detector offset: the detector is shifted sideways and geo.offDetector is not set.
    if wanted("offset"):
        g = G()
        g.offDetector = np.array([0.0, 6.0])  # (v, u) in Python: 6 mm along u
        p = tigre.Ax(head, g, angles)
        figure(out, "offset", [
            ("Correct", axial(good), None),
            ("Detector 6 mm off, not in geo", axial(algs.fdk(p, G(), angles)), None),
            ("Same data, offDetector set", axial(algs.fdk(p, g, angles)), None),
        ])

    # 3. Angles: the scan turned the other way, or started elsewhere.
    if wanted("angles"):
        figure(out, "angles", [
            ("Correct", axial(good), None),
            ("Angles in the wrong direction", axial(algs.fdk(proj, G(), -angles)), None),
            ("Angles start 30° off", axial(algs.fdk(proj, G(), angles + np.pi / 6)), None),
        ])

    # 4. Log transform and flat field: raw intensities instead of line integrals.
    if wanted("log"):
        i0 = 60000.0
        intensity = (i0 * np.exp(-proj)).astype(np.float32)
        no_log = algs.fdk(intensity / i0, G(), angles)
        no_flat = algs.fdk(-np.log(intensity), G(), angles)  # log, but no division by I0
        figure(out, "log", [
            ("Correct", axial(good), None),
            ("No log (I/I0 given)", axial(no_log), auto_window(axial(no_log))),
            ("Log without flat field", axial(no_flat), auto_window(axial(no_flat))),
        ])

    # 5. Short scan: 180° plus the fan angle, with and without Parker weights.
    if wanted("shortscan"):
        fan = 2 * np.arctan(geo.sDetector[1] / 2 / geo.DSD)
        short = np.linspace(0, np.pi + fan, 240, endpoint=False)
        p = tigre.Ax(head, G(), short)
        no_parker = algs.fdk(p, G(), short, parker=False)
        diff = axial(no_parker) - axial(good)
        figure(out, "shortscan", [
            ("Short scan, no Parker weights", axial(no_parker), None),
            ("Its difference from a full scan", diff, diff_window(diff)),
            ("Short scan, Parker weights", axial(algs.fdk(p, G(), short, parker=True)), None),
        ])

    # 6. Too few projections: streaks in FDK; an iterative method copes better.
    if wanted("fewangles"):
        few = np.linspace(0, 2 * np.pi, 40, endpoint=False)
        p = tigre.Ax(head, G(), few)
        figure(out, "fewangles", [
            ("360 projections, FDK", axial(good), None),
            ("40 projections, FDK", axial(algs.fdk(p, G(), few)), None),
            ("40 projections, OS-ASD-POCS", axial(algs.os_asd_pocs(p, G(), few, 100)), None),
        ])

    # 7. Iterations on noisy data: too few blur, too many amplify the noise.
    if wanted("iterations"):
        noisy = CTnoise.add(proj, Poisson=2e4, Gaussian=np.array([0, 10])).astype(np.float32)
        figure(out, "iterations", [
            ("CGLS, 3 iterations", axial(algs.cgls(noisy, G(), angles, 3)), None),
            ("CGLS, 20 iterations", axial(algs.cgls(noisy, G(), angles, 20)), None),
            ("CGLS, 150 iterations", axial(algs.cgls(noisy, G(), angles, 150)), None),
        ])

    # 8. Truncation: the object is wider than the detector sees.
    if wanted("truncation"):
        g = G()
        g.nDetector = np.array([N, 150])
        g.sDetector = g.nDetector * g.dDetector
        p = tigre.Ax(head, g, angles)
        figure(out, "truncation", [
            ("Whole object on the detector", axial(good), None),
            ("Detector too narrow, FDK", axial(algs.fdk(p, g, angles)), None),
        ])

    # 9. Distances: DSO or DSD wrong by 5 %.
    if wanted("distances"):
        g = G()
        g.DSO = geo.DSO * 1.10
        wrong = axial(algs.fdk(proj, g, angles))
        diff = wrong - axial(good)
        figure(out, "distances", [
            ("Correct", axial(good), None),
            ("DSO 10 % too large", wrong, None),
            ("Difference", diff, diff_window(diff)),
        ])

    # 10. Detector roll: the detector is rotated in its own plane and rotDetector is not set.
    if wanted("roll"):
        g = G()
        g.rotDetector = np.array([np.deg2rad(3.0), 0.0, 0.0])
        p = tigre.Ax(head, g, angles)
        figure(out, "roll", [
            ("Correct (coronal)", coronal(good), None),
            ("Detector rolled 3°, not in geo", coronal(algs.fdk(p, G(), angles)), None),
            ("Same data, rotDetector set", coronal(algs.fdk(p, g, angles)), None),
        ])

    # 11. Detector axes swapped when loading: (v, u) read as (u, v).
    if wanted("axes"):
        swapped = np.ascontiguousarray(np.transpose(proj, (0, 2, 1)))
        bad = algs.fdk(swapped, G(), angles)
        figure(out, "axes", [
            ("Correct", axial(good), None),
            ("Projections transposed (u and v swapped)", axial(bad), auto_window(axial(bad))),
        ])

    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
