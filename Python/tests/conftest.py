"""Matplotlib compatibility helpers shared by TIGRE's tests."""
import numpy as np


def poly_vertices(collection):
    """Return a Poly3DCollection's vertices as a (coordinates, N) array.

    Matplotlib 3.10 and earlier keep them in ``_vec``; 3.11 keeps the faces in ``_faces``.
    """
    if hasattr(collection, "_vec"):
        return collection._vec
    return np.asarray(collection._faces).reshape(-1, 3).T
