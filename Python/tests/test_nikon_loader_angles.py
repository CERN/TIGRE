"""NikonDataLoader reads the angles from each kind of file a Nikon scan can have.

The ``.ang`` file of some scanners (for example the SophiaBeads datasets,
https://doi.org/10.5281/zenodo.16474) has tab-separated lines such as
``1:\t005.716``; TIGRE 3.3.1 and earlier passed the text to ``math.radians``
and failed with "must be real number, not str".
"""
import numpy as np
import pytest

from tigre.utilities.io.NikonDataLoader import readXtekctGeometry

XTEKCT = """[XTekCT]
Name=sample
VoxelsX=8
VoxelsY=8
VoxelsZ=4
VoxelSizeX=0.1
VoxelSizeY=0.1
VoxelSizeZ=0.1
SrcToObject=80.0
SrcToDetector=1000.0
DetectorPixelsX=10
DetectorPixelsY=6
DetectorPixelSizeX=0.2
DetectorPixelSizeY=0.2
DetectorOffsetX=0.0
DetectorOffsetY=0
CentreOfRotationTop=0.0
WhiteLevel=60000.0
Projections=3
InitialAngle=0
AngularStep=5.0
"""


def _scan(tmp_path, name=None, text=None):
    (tmp_path / "sample.xtekct").write_text(XTEKCT)
    if name:
        # newline="" keeps the Windows line endings the scanner writes
        with open(tmp_path / name, "w", newline="") as f:
            f.write(text)
    return str(tmp_path)


@pytest.mark.parametrize("line_end", ["\n", "\r\n"])
def test_ang_file_with_tabs(tmp_path, line_end):
    text = line_end.join(["Proj\tAngle(deg)", "1:\t005.716", "2:\t000.000", "3:\t-005.716"]) + line_end
    _, _, angles = readXtekctGeometry(_scan(tmp_path, "sample.ang", text))
    np.testing.assert_allclose(angles, np.radians([5.716, 0.0, -5.716]))


def test_ang_file_with_colons(tmp_path):
    text = "Angles\n1: 0.0\n2: 5.0\n3: 10.0\n"
    _, _, angles = readXtekctGeometry(_scan(tmp_path, "sample.ang", text))
    np.testing.assert_allclose(angles, np.radians([0.0, 5.0, 10.0]))


def test_ctdata_file(tmp_path):
    text = "header\nheader\nheader\n1\t0.0\n2\t5.0\n3\t10.0\n"
    _, _, angles = readXtekctGeometry(_scan(tmp_path, "sample_ctdata.txt", text))
    np.testing.assert_allclose(angles, -np.radians([0.0, 5.0, 10.0]))


def test_no_angle_file(tmp_path):
    _, geo, angles = readXtekctGeometry(_scan(tmp_path))
    np.testing.assert_allclose(angles, -np.radians([0.0, 5.0, 10.0]))
    assert geo.DSO == 80.0 and geo.DSD == 1000.0
