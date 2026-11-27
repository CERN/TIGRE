from pathlib import Path
import warnings

import numpy as np

from tigre.utilities.geometry import Geometry

# Port of MATLAB/Utilities/IO/Dicom/dicomDataLoader.m, a modification of code written for:
# Hatamikia, S., Biguri, A., Kronreif, G., et al. (2020), Optimization for customized
# trajectories in cone beam computed tomography. Med. Phys., 47: 4786-4799.
# https://doi.org/10.1002/mp.14403 . Please cite it, especially if you are using C-arm CBCT.

_REQUIRED = (
    "DistanceSourceToDetector",
    "DistanceSourceToPatient",
    "Rows",
    "Columns",
    "PositionerPrimaryAngle",
    "PositionerSecondaryAngle",
    "ImagerPixelSpacing",
)
_INCREMENTS = ("PositionerPrimaryAngleIncrement", "PositionerSecondaryAngleIncrement")
_SHUTTER = (
    "ShutterLeftVerticalEdge",
    "ShutterRightVerticalEdge",
    "ShutterUpperHorizontalEdge",
    "ShutterLowerHorizontalEdge",
)
_CONTACT = (
    "If you have data not supported by this loader that you want supported, contact "
    "tigre.toolbox@gmail.com and we can work together in improving this code to load your data."
)


def DicomDataLoader(filepath):
    """Load a DICOM C-arm CBCT dataset: projections, geometry and angles.

    Reads every ``.dcm`` file in a folder. A standard scan is one multi-frame
    file, with the angle of each frame in
    ``PositionerPrimaryAngleIncrement`` / ``PositionerSecondaryAngleIncrement``;
    a scan made of one-shot acquisitions is one file per projection. The
    projections are converted to line integrals, ready to reconstruct. Tested
    on Philips C-arm (Allura) data only. Needs the ``pydicom`` package
    (``pip install pydicom``). See :doc:`/user-guide/loading-data`.

    Parameters
    ----------
    filepath : str
        Folder with the ``.dcm`` files.

    Returns
    -------
    proj : numpy.ndarray
        Projections (line integrals), float32, shape ``(N, v, u)``.
    geo : tigre.utilities.geometry.Geometry
        The scanner geometry, with an image size estimated from the detector.
        ``None`` if the headers lack the information to build it; ``proj`` is
        then the raw pixel data.
    angles : numpy.ndarray
        Projection angles in radians: shape ``(N,)``, or ``(N, 3)`` if the
        C-arm also moved in its secondary (cranial/caudal) angle. ``None`` if
        ``geo`` is.
    headers : list of pydicom.dataset.FileDataset
        The DICOM headers of the files used, without pixel data.
    """
    try:
        import pydicom
    except ImportError as err:
        raise ImportError("DicomDataLoader needs pydicom: pip install pydicom") from err

    folder = Path(filepath)
    files = sorted(folder.glob("*.dcm"))
    if not files:
        raise ValueError(f"No .dcm files found in folder: {folder}")
    headers = [pydicom.dcmread(f, stop_before_pixels=True) for f in files]

    last = headers[-1]
    if "Philips" not in str(getattr(last, "Manufacturer", "")):
        warnings.warn(
            "This loader has only been tested for Philips scanners, which your DICOM data does not "
            "seem to be. " + _CONTACT + " Attempting to read the data anyway."
        )
    if str(getattr(last, "PositionerType", "")).upper() != "CARM":
        warnings.warn(
            "This loader has only been tested with C-arm (CARM) scanners, which your DICOM data does "
            "not seem to be. " + _CONTACT + " Attempting to read the data anyway."
        )

    off_u = off_v = 0.0
    if len(files) > 1:
        # One projection per file. Some headers may be corrupted: drop those.
        headers, files = _check_headers(headers, files)
        if not headers:
            raise ValueError("No DICOM file has the fields needed to build a geometry.")
        proj = np.stack([_first_frame(pydicom.dcmread(f)) for f in files]).astype(np.float32)
        DSD = np.array([float(h.DistanceSourceToDetector) for h in headers])
        DSO = np.array([float(h.DistanceSourceToPatient) for h in headers])
        primary = np.array([float(h.PositionerPrimaryAngle) for h in headers])  # RAO/LAO
        secondary = np.array([float(h.PositionerSecondaryAngle) for h in headers])  # CRA/CAU
    else:
        # A standard scan: one multi-frame file.
        valid, _ = _check_headers(headers, files, _INCREMENTS)
        data = pydicom.dcmread(files[0])
        if not valid:
            warnings.warn(
                "DICOM header lacks required information to construct a geometry, "
                "returning projection data only."
            )
            return data.pixel_array, None, None, headers
        h = valid[0]
        proj = _all_frames(data).astype(np.float32)
        DSD = np.array([float(h.DistanceSourceToDetector)])
        DSO = np.array([float(h.DistanceSourceToPatient)])
        primary = float(h.PositionerPrimaryAngle) + _as_array(h.PositionerPrimaryAngleIncrement)
        secondary = float(h.PositionerSecondaryAngle) + _as_array(h.PositionerSecondaryAngleIncrement)
        if all(hasattr(h, f) for f in _SHUTTER):
            # Keep only the pixels inside the shutter (DICOM edges are 1-based and inclusive).
            rows, cols = proj.shape[1], proj.shape[2]
            top, bottom = int(h.ShutterUpperHorizontalEdge), int(h.ShutterLowerHorizontalEdge)
            left, right = int(h.ShutterLeftVerticalEdge), int(h.ShutterRightVerticalEdge)
            proj = proj[:, top - 1 : bottom, left - 1 : right]
            # The cropped detector is no longer centred: move it back to where its pixels were.
            # u grows with the column and v with the row (TIGRE's row 0 is the -z end of the detector).
            d_v, d_u = (float(x) for x in h.ImagerPixelSpacing)
            off_u = ((left + right) / 2 - (cols + 1) / 2) * d_u
            off_v = ((top + bottom) / 2 - (rows + 1) / 2) * d_v

    if np.any(proj == 0):
        warnings.warn(
            "Zeroes detected in projection data, adding +1 to apply the Beer-Lambert law. This is "
            "caused either by very high attenuation, to the point of blocking x-rays totally, or by "
            "the detector being cropped digitally (horizontal/vertical zero bands). In the first case, "
            "you may consider scanning your data again; in the second, cropping the projections and "
            "changing the geo.*Detector fields. It can also be a single rogue pixel, and you may "
            "continue safely."
        )
        proj = proj + 1
    proj = -np.log(proj / (proj.max() + 1)).astype(np.float32)

    # Geometry. ImagerPixelSpacing is (row, column) spacing, i.e. (v, u): TIGRE's Python order.
    geo = Geometry()
    geo.mode = "cone"
    geo.accuracy = 0.5
    geo.dDetector = np.array([float(x) for x in headers[0].ImagerPixelSpacing])
    geo.nDetector = np.array(proj.shape[1:])
    geo.sDetector = geo.dDetector * geo.nDetector
    geo.offDetector = np.array([off_v, off_u])
    geo.DSD = DSD[0] if np.unique(DSD).size == 1 else DSD
    geo.DSO = DSO[0] if np.unique(DSO).size == 1 else DSO  # very unlikely in a medical scan

    # Estimate an acceptable image size, as in the MATLAB loader.
    d_v, d_u = geo.dDetector
    mag = np.mean(DSO) / np.mean(DSD)
    geo.dVoxel = np.array([d_v, d_u, d_u]) * mag
    n_xy = np.ceil(geo.nDetector[1] + abs(geo.offDetector[1]) / d_u)
    geo.nVoxel = np.array([geo.nDetector[0], n_xy, n_xy]).astype(int)
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.offOrigin = np.array([0, 0, 0])
    geo.rotDetector = np.array([0, 0, 0])

    if np.all(secondary == 0):
        angles = -primary * np.pi / 180
    else:
        angles = np.stack([-primary, -secondary, np.zeros_like(primary)], axis=1) * np.pi / 180

    print("Projection data loaded.")
    return proj, geo, angles, headers


def _check_headers(headers, files, extra=()):
    """Drop the headers (and files) that lack a field needed to build the geometry."""
    required = _REQUIRED + tuple(extra)
    keep = [i for i, h in enumerate(headers) if all(hasattr(h, f) for f in required)]
    if len(keep) < len(headers):
        warnings.warn(
            f"Some dicom files did not contain required fields: {len(headers) - len(keep)} files ignored."
        )
    return [headers[i] for i in keep], [files[i] for i in keep]


def _as_array(value):
    """A DICOM value that may hold one number or several, as a float array."""
    if isinstance(value, (list, tuple)) or hasattr(value, "__iter__") and not isinstance(value, str):
        return np.array([float(v) for v in value])
    return np.array([float(value)])


def _all_frames(data):
    """Pixel data as (frames, rows, columns), keeping the first sample of colour data."""
    arr = data.pixel_array
    frames = int(getattr(data, "NumberOfFrames", 1) or 1)
    samples = int(getattr(data, "SamplesPerPixel", 1) or 1)
    if samples > 1:
        arr = arr[..., 0]
    if frames == 1:
        arr = arr[np.newaxis]
    return arr


def _first_frame(data):
    """The first frame of a file's pixel data, as (rows, columns)."""
    return _all_frames(data)[0]
