import importlib.util
import struct
from pathlib import Path

import numpy as np

XIM_PATH = Path(__file__).resolve().parents[1] / "tigre" / "utilities" / "io" / "varian" / "xim.py"


def _load_xim():
    spec = importlib.util.spec_from_file_location("xim", XIM_PATH)
    xim = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(xim)
    return xim.XIM


def test_uncompressed_xim_exposes_pixel_array(tmp_path):
    xim_path = tmp_path / "uncompressed.xim"
    pixels = np.array([[-1, 0], [1, 1024]], dtype="<i4")
    xim_path.write_bytes(
        b"XIMFILE\x00"
        + struct.pack("<6i", 1, 2, 2, 32, 4, 0)
        + struct.pack("<i", pixels.nbytes)
        + pixels.tobytes()
        + struct.pack("<2i", 0, 0)
    )

    image = _load_xim()(xim_path)

    assert image.array.shape == (2, 2)
    np.testing.assert_array_equal(image.array, pixels)
