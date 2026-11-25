"""Make the figures and numbers of the website's real-data tutorial (SophiaBeads).

The SophiaBeads datasets are on Zenodo (https://doi.org/10.5281/zenodo.16474). With
TIGRE installed and a GPU, this downloads the smallest one (64 projections, a
0.47 GB zip, 0.5 GB unzipped) into ``--data`` the first time, unzips it, and runs
the tutorial. It needs about 1 GB of free disk space while unzipping and about
6 GB of RAM. The data are not part of TIGRE; please cite them if you use them
(see the Zenodo record).

    python make_sophiabeads_tutorial.py --out <output folder>

To use a copy you already have, give its folder instead:

    python make_sophiabeads_tutorial.py <dataset folder> --out <output folder>

(``zenodo_get 16474 -g "*64*"`` from the ``zenodo-get`` package downloads the same
file by hand.)

It follows the tutorial's steps (load, look, check the geometry on a few slices,
centre of rotation, full reconstruction with FDK and CGLS) and
writes, to the output folder, one PNG per figure and ``numbers.json`` with the
values the page quotes (geometry, timings, the centre of rotation).
"""

import argparse
import copy
import json
import pathlib
import time
import urllib.request
import zipfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import tigre  # noqa: E402
import tigre.algorithms as algs  # noqa: E402
from tigre.utilities.io import NikonDataLoader  # noqa: E402


ZENODO_RECORD = "16474"  # SophiaBeads Datasets


def download_sophiabeads(data_dir, n_projections=64):
    """Download and unzip the SophiaBeads dataset with n_projections projections.

    Returns the folder that holds its .xtekct file. Does nothing if it is
    already there.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    found = sorted(data_dir.rglob("*.xtekct"))
    found = [f for f in found if f"_{n_projections}" in f.parent.name] or found
    if found:
        return found[0].parent
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{ZENODO_RECORD}") as response:
        files = json.load(response)["files"]
    wanted = [f for f in files if f"_{n_projections}" in f["key"] and f["key"].endswith(".zip")]
    if not wanted:
        raise RuntimeError("No file with %d projections in Zenodo record %s: %s"
                           % (n_projections, ZENODO_RECORD, [f["key"] for f in files]))
    entry = wanted[0]
    target = data_dir / entry["key"]
    tmp = target.with_suffix(".part")
    size = entry.get("size", 0)
    if not target.exists() and not (tmp.exists() and tmp.stat().st_size == size):
        print(f"Downloading {entry['key']} ({size / 1e9:.2f} GB) from Zenodo into {data_dir}.\n"
              f"You need about {2 * size / 1e9 + 0.1:.0f} GB of free disk space while it is unzipped, "
              "and about 6 GB of RAM for the tutorial.")
        with urllib.request.urlopen(entry["links"]["self"]) as response, open(tmp, "wb") as out:
            done = 0
            while chunk := response.read(1 << 22):
                out.write(chunk)
                done += len(chunk)
                print(f"\r  {done / 1e9:.2f} GB", end="", flush=True)
        print()
    if not target.exists():
        # On Windows the antivirus may still hold the new file for a moment.
        for attempt in range(20):
            try:
                tmp.replace(target)
                break
            except PermissionError:
                if attempt == 19:
                    raise
                time.sleep(1)
    print(f"Unzipping {target.name}...")
    with zipfile.ZipFile(target) as archive:
        archive.extractall(data_dir)
    try:
        target.unlink()  # the unzipped copy is all the tutorial needs
    except PermissionError:
        print(f"Could not delete {target}; you can delete it yourself.")
    found = sorted(data_dir.rglob("*.xtekct"))
    if not found:
        raise RuntimeError(f"No .xtekct file in {target.name}")
    return found[0].parent


def load_nikon(folder):
    """NikonDataLoader, with a work-around for TIGRE 3.3.1 and earlier.

    Those versions read the angles of a .ang file (as in SophiaBeads) as text
    and fail with "must be real number, not str".
    """
    try:
        return NikonDataLoader(folder)
    except TypeError as error:
        if "real number" not in str(error):
            raise
        import math
        import types

        module = __import__(NikonDataLoader.__module__, fromlist=["math"])
        module.math = types.SimpleNamespace(radians=lambda x: math.radians(float(x)))
        try:
            return NikonDataLoader(folder)
        finally:
            module.math = math


def as_list(x):
    return np.asarray(x, dtype=float).round(4).tolist()


def save(fig, out, name):
    path = out / f"{name}.png"
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def show_slices(out, name, panels, window):
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels), 3.9))
    for ax, (title, img) in zip(np.atleast_1d(axes), panels):
        ax.imshow(img, cmap="gray", vmin=window[0], vmax=window[1])
        ax.set_title(title, fontsize=11)
        ax.axis("off")
    fig.tight_layout()
    save(fig, out, name)


def slab(geo, n_slices):
    """The same geometry, but only n_slices slices around the centre."""
    g = copy.deepcopy(geo)
    g.nVoxel = np.array([n_slices, geo.nVoxel[1], geo.nVoxel[2]])
    g.sVoxel = g.nVoxel * g.dVoxel
    return g


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dataset", type=pathlib.Path, nargs="?",
                        help="unzipped SophiaBeads folder (with the .xtekct file); "
                             "default: download it into --data")
    parser.add_argument("--data", type=pathlib.Path, default=pathlib.Path("sophiabeads_data"),
                        help="where to download the dataset (default: ./sophiabeads_data)")
    parser.add_argument("--out", type=pathlib.Path, default=pathlib.Path("sophiabeads_figures"))
    args = parser.parse_args()
    dataset = args.dataset or download_sophiabeads(args.data)
    print(f"Dataset: {dataset}")
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    numbers = {}

    # 1. Load: projections as line integrals, the geometry and the angles.
    t = time.time()
    proj, geo, angles = load_nikon(str(dataset))
    numbers["load_seconds"] = round(time.time() - t, 1)
    numbers["proj_shape"] = list(proj.shape)
    numbers["proj_range"] = [float(np.percentile(proj, 0.1)), float(np.percentile(proj, 99.9))]
    numbers["n_angles"] = len(angles)
    numbers["angles_deg_first_last"] = [float(np.rad2deg(angles[0])), float(np.rad2deg(angles[-1]))]
    numbers["geo"] = {k: as_list(getattr(geo, k)) for k in
                      ("DSO", "DSD", "nDetector", "dDetector", "offDetector", "nVoxel", "dVoxel", "sVoxel")
                      if hasattr(geo, k)}
    numbers["geo"]["COR_from_file"] = as_list(getattr(geo, "COR", 0.0))
    print(json.dumps(numbers, indent=1))

    # 2. Look at one projection.
    fig, ax = plt.subplots(figsize=(4.5, 4.2))
    im = ax.imshow(proj[0], cmap="gray")
    ax.set_title("Projection 1 (line integrals)", fontsize=11)
    ax.axis("off")
    fig.colorbar(im, ax=ax, fraction=0.046)
    save(fig, out, "projection")

    # 3. FDK on a few slices around the centre, as the file's geometry says.
    test = slab(geo, 8)
    t = time.time()
    img_file = algs.fdk(proj, copy.deepcopy(test), angles)
    numbers["fdk_slab_seconds"] = round(time.time() - t, 1)
    mid = img_file.shape[0] // 2
    lo, hi = np.percentile(img_file[mid], [1, 99.8])
    window = (float(lo), float(hi))
    numbers["window"] = window

    # 4. Centre of rotation, estimated from the data.
    cor = float(tigre.computeCOR(proj, copy.deepcopy(geo), angles))
    numbers["COR_estimated"] = round(cor, 4)
    g = copy.deepcopy(test)
    g.COR = cor
    img_cor = algs.fdk(proj, g, angles)
    crop = slice(img_cor.shape[1] // 3, 2 * img_cor.shape[1] // 3)
    show_slices(out, "cor", [
        ("Geometry from the file", img_file[mid]),
        (f"geo.COR = {cor:.3f} mm", img_cor[mid]),
    ], window)
    show_slices(out, "cor_zoom", [
        ("Geometry from the file (zoom)", img_file[mid][crop, crop]),
        (f"geo.COR = {cor:.3f} mm (zoom)", img_cor[mid][crop, crop]),
    ], window)
    geo.COR = cor

    # 5. The whole object: FDK and CGLS. The file's grid is
    #    2000^3 voxels (32 GB as float32); a quarter of the resolution (500^3,
    #    voxels 4x larger, same field of view) fits most computers and GPUs.
    half = copy.deepcopy(geo)
    half.nVoxel = (geo.nVoxel // 4).astype(int)
    half.dVoxel = geo.sVoxel / half.nVoxel
    half.sVoxel = half.nVoxel * half.dVoxel
    numbers["full_nVoxel"] = as_list(half.nVoxel)
    results = {}
    for label, run in [
        ("FDK", lambda: algs.fdk(proj, copy.deepcopy(half), angles)),
        ("CGLS, 15 iterations", lambda: algs.cgls(proj, copy.deepcopy(half), angles, 15)),
    ]:
        t = time.time()
        results[label] = run()
        numbers[f"seconds {label}"] = round(time.time() - t, 1)
        print(f"{label}: {numbers[f'seconds {label}']} s")
    zmid = half.nVoxel[0] // 2
    show_slices(out, "full_axial", [(k, v[zmid]) for k, v in results.items()], window)
    show_slices(out, "full_coronal", [(k, v[:, v.shape[1] // 2, :]) for k, v in results.items()], window)

    # 6. Save the reconstruction (NumPy; any format works, e.g. tifffile.imwrite).
    np.save(out / "sophiabeads_cgls.npy", results["CGLS, 15 iterations"])

    (out / "numbers.json").write_text(json.dumps(numbers, indent=1))
    print(f"wrote {out / 'numbers.json'}")


if __name__ == "__main__":
    main()
