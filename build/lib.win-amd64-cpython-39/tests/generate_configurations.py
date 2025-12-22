import os

import numpy as np
import tigre

nVoxel = np.array([256, 256, 256])  # noqa: N816
niter = 20
dirname = os.path.dirname(__file__)
keywords = dict(blocksize=20, hyper=9.391252e06)


def save_config(dict, name):
    for kw in dict:
        if kw not in ["nproj", "geo", "angles", "niter", "kwargs"]:
            raise KeyError()
    np.save(os.path.join(dirname, name), dict)


if __name__ == "__main__":

    geo = tigre.geometry(mode="cone", nVoxel=nVoxel, default=True)
    # (V,U) number of pixels        (px)
    geo.nDetector = np.array((256, 256))
    # size of each pixel            (mm)
    geo.dDetector = np.array((0.8, 0.8)) * 2
    geo.sDetector = geo.nDetector * geo.dDetector

    save_config(
        dict(nproj=100, niter=20, geo=geo, angles=np.linspace(0, 2 * np.pi, 100), kwargs=keywords),
        "configuration1.npy",
    )

    save_config(
        dict(
            nproj=100,
            niter=20,
            geo=tigre.geometry(mode="parallel", default=True, nVoxel=np.array([256, 256, 256])),
            angles=np.linspace(0, 2 * np.pi, 100),
            kwargs=keywords,
        ),
        "configuration2.npy",
    )
    geo = tigre.geometry(mode="cone", nVoxel=np.array([240, 248, 256]), default=True)
    geo.nDetector = np.array((248, 256))
    geo.dDetector = np.array((0.8, 0.8)) * 2
    geo.sDetector = geo.nDetector * geo.dDetector
    save_config(
        dict(nproj=100, niter=20, geo=geo, angles=np.linspace(0, 2 * np.pi, 100), kwargs=keywords),
        "configuration3.npy",
    )

    save_config(
        dict(
            nproj=100,
            niter=20,
            geo=tigre.geometry(mode="parallel", default=True, nVoxel=np.array([240, 248, 256])),
            angles=np.linspace(0, 2 * np.pi, 100),
            kwargs=keywords,
        ),
        "configuration4.npy",
    )
