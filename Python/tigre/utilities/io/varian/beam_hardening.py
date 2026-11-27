"""Beam-hardening correction for Varian TrueBeam CBCT, ported from MATLAB.

Port of MATLAB/Utilities/IO/VarianCBCT/BHCorrection.m and the BH_*.m files it
calls (Yi Du, 2021-2022). The correction models the polychromatic spectrum
after the kV filter and the bowtie, builds the log signal that water of each
thickness would give behind each bowtie path length, and remaps every pixel
of the log-normalised projections to the linear (monochromatic) signal of the
equivalent water thickness.

Differences from MATLAB, none of which change the model:
- the path length through the bowtie is computed by tracing each detector ray
  through the bowtie profile in NumPy, instead of voxelising the bowtie
  (4002 x 4002 x 102) and forward projecting it with ``Ax``;
- the bowtie profile is fitted with SciPy's smoothing spline instead of the
  Curve Fitting Toolbox's (the smoothing parameter is chosen differently);
- the ``a*x`` fits are solved in closed form (least squares);
- negative and unmappable pixels are filled by a least-squares smooth
  interpolation (second differences), as ``inpaint_nans`` does.
"""

from __future__ import annotations

import glob
import os
import xml.etree.ElementTree as ET

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import make_smoothing_spline
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from tqdm import tqdm

from tigre.utilities.geometry import Geometry
from tigre.utilities.io.varian.utils import XML, PathLike

# Height of the bowtie (TrueBeam Technical Reference), mm: rays above or below it miss it.
BOWTIE_HEIGHT = 55.0
# Thickness of the CsI scintillator, mm. Cancels out in the log signal.
SCINTILLATOR_THICKNESS = 0.6


def _read_bh_xml(filepath: PathLike) -> ET.Element:
    """Reads Calibrations/ASC/Factory/Calibration.xml from the scan folder."""
    file = os.path.join(filepath, "Calibrations", "ASC", "Factory", "Calibration.xml")
    if not os.path.isfile(file):
        files = glob.glob(os.path.join(filepath, "Calibrations", "ASC*", "Factory", "Calibration.xml"))
        if len(files) != 1:
            raise RuntimeError("Beam hardening calibration file (Calibrations/ASC/Factory/Calibration.xml) not found.")
        file = files[0]
    return ET.parse(file).getroot()


class BHParams(XML):
    """Beam-hardening calibration (spectrum, filter, bowtie, materials) for one scan."""

    def __init__(self, filepath: PathLike, voltage: str, kv_filter: str, bowtie: str, xml_reader=_read_bh_xml) -> None:
        super().__init__(filepath, xml_reader)
        res = "CalibrationResults/"
        self.energy_bin = float(self._get_field(res + "EnergyBinWidth").text)

        # Source spectrum at the scan voltage
        for spec in self.root.findall(res + "Spectra/SpectrumProperties", self.ns):
            if self._text(spec, "Voltage").lower() == voltage.strip().lower():
                self.voltage = float(voltage)
                self.flux = self._floats(spec, "Flux")
                break
        else:
            raise RuntimeError(f"No spectrum for {voltage} kV in the beam hardening calibration.")
        self.source_spec = self.flux / self.flux.max()
        self.kV = np.arange(1, self.voltage + self.energy_bin / 2, self.energy_bin)

        # kV filter
        self.filter_name = kv_filter
        self.filter_material = None
        self.filter_thickness = 0.0
        for f in self.root.findall(res + "Filters/FilterProperties", self.ns):
            if self._text(f, "Id").lower() == kv_filter.strip().lower():
                self.filter_material = self._text(f, "Material")
                self.filter_thickness = float(self._text(f, "Thickness"))
                break

        # Bowtie: distance from the source (mm), material, profile (thickness, mm, vs u, mm)
        self.bowtie_name = bowtie
        self.bowtie_material = None
        for b in self.root.findall(res + "Bowties/BowtieProperties", self.ns):
            if b.find("Id", self.ns) is None:  # some items are comments
                continue
            if self._text(b, "Id").lower() == bowtie.strip().lower():
                self.bowtie_distance = float(self._text(b, "SourceDist"))
                if b.find("Coordinate", self.ns) is not None:
                    self.bowtie_material = self._text(b, "Material")
                    self.bowtie_u = self._floats(b, "Coordinate")
                    self.bowtie_thickness = self._floats(b, "Profile")
        if self.bowtie_material is None:
            raise RuntimeError(f"No profile for bowtie '{bowtie}' in the beam hardening calibration.")

        # Calibration object (water by default)
        self.object_material = self._get_field(res + "ObjectMaterialId").text
        self.object_thickness_max = float(self._get_field(res + "ObjectThicknessMax").text)

        # Attenuation coefficients (/mm) of the materials, per energy bin
        self.filter_ac = self.object_ac = self.bowtie_ac = self.scintillator_ac = None
        for m in self.root.findall(res + "Materials/MaterialProperties", self.ns):
            name = self._text(m, "Id")
            if self.filter_material is not None and name.lower() == self.filter_material.lower():
                self.filter_ac = self._floats(m, "AttenuationCoefficient")
            if name.lower() == self.object_material.lower():
                self.object_ac = self._floats(m, "AttenuationCoefficient")
            if name.lower() == self.bowtie_material.lower():
                self.bowtie_ac = self._floats(m, "AttenuationCoefficient")
            if "cesium" in name.lower():  # scintillator (PaxScan 4030CB: CsI:Tl)
                self.scintillator_ac = self._floats(m, "EnergyAbsorptionCoefficient")
        for attr, what in [("object_ac", self.object_material), ("bowtie_ac", self.bowtie_material),
                           ("scintillator_ac", "cesium (scintillator)")]:
            if getattr(self, attr) is None:
                raise RuntimeError(f"No attenuation coefficients for {what} in the beam hardening calibration.")

    def _text(self, elem: ET.Element, field: str) -> str:
        sub = elem.find(field, self.ns)
        if sub is None or sub.text is None:
            raise RuntimeError(f"{field} not found in {self.name}")
        return sub.text.strip()

    def _floats(self, elem: ET.Element, field: str) -> NDArray:
        sub = elem.find(field, self.ns)
        if sub is None:
            raise RuntimeError(f"{field} not found in {self.name}")
        return np.array([float(f.text) for f in sub.findall("float", self.ns)])

    def filtered_spectrum(self) -> NDArray:
        """Source spectrum after the kV filter, normalised to 1 (BH_SpectrumFilter.m)."""
        if "none" in self.filter_name.lower() or self.filter_ac is None:
            return self.source_spec
        spec = self.source_spec * np.exp(-self.filter_thickness * self.filter_ac[: len(self.source_spec)])
        return spec / spec.max()


def bowtie_path_lengths(geometry: Geometry, bh: BHParams, n_steps: int = 4000) -> NDArray:
    """Length (mm) of each detector ray inside the bowtie, shape (nv, nu).

    The bowtie is the region between ``bowtie_distance`` and
    ``bowtie_distance + thickness(y)`` from the source, where ``y`` is the
    lateral (u) coordinate at the bowtie and ``thickness`` its fitted profile,
    and less than ``BOWTIE_HEIGHT / 2`` from the central plane: the body that
    BH_SpectrumBowtieLUT.m voxelises and forward projects.
    """
    DSD = float(np.mean(geometry.DSD))
    SBD = bh.bowtie_distance
    d_v, d_u = np.asarray(geometry.dDetector, dtype=np.float64)[:2]
    n_v, n_u = (int(n) for n in np.asarray(geometry.nDetector)[:2])
    off = np.asarray(geometry.offDetector, dtype=np.float64).reshape(-1, 2)[0] if geometry.offDetector is not None else np.zeros(2)
    u = (np.arange(n_u) - n_u / 2 + 0.5) * d_u + off[1]
    v = (np.arange(n_v) - n_v / 2 + 0.5) * d_v + off[0]

    # Bowtie profile: smoothing spline through the calibration points
    order = np.argsort(bh.bowtie_u)
    bu, bt = bh.bowtie_u[order], bh.bowtie_thickness[order]
    bu, keep = np.unique(bu, return_index=True)
    bt = bt[keep]
    profile = make_smoothing_spline(bu, bt)
    ss_res = np.sum((profile(bu) - bt) ** 2)
    ss_tot = np.sum((bt - bt.mean()) ** 2)
    print(f"Bowtie fitting goodness (R^2): {1 - ss_res / ss_tot if ss_tot > 0 else 1.0:.6f}")
    t_max = np.ceil(bt.max())

    # March along each ray from the source side of the bowtie to t_max behind it.
    depth = (np.arange(n_steps) + 0.5) * (t_max / n_steps)  # mm past the bowtie front face
    dist = SBD + depth  # distance from the source, along the central ray
    y = np.outer(dist, u) / DSD  # lateral position at that distance, (n_steps, n_u)
    thick = np.clip(profile(y), 0, t_max)
    inside = (depth[:, None] < thick).astype(np.float64)
    # Rays also leave the bowtie above/below it: keep only depths with |z| <= height/2.
    cum = np.vstack([np.zeros((1, n_u)), np.cumsum(inside, axis=0)]) * (t_max / n_steps)
    with np.errstate(divide="ignore"):
        dist_max = np.where(v == 0, np.inf, (BOWTIE_HEIGHT / 2) * DSD / np.abs(v))
    k = np.searchsorted(dist, dist_max, side="right")  # depth samples inside the height
    length = cum[k]  # (n_v, n_u), measured along the beam axis
    # ... and along the ray itself
    length *= np.sqrt(DSD**2 + u[None, :] ** 2 + v[:, None] ** 2) / DSD
    return length


def water_lookup(bh: BHParams, bowtie_lengths: NDArray, n_bowtie: int = 2000):
    """Log signal of each water thickness behind each bowtie length (BH_SpectrumBowtieLUT.m,
    BH_ObjectCalibLUT.m) and the slope of its linear part (BH_RemappingFunc.m).

    Returns:
        bowtie_sl (n_bowtie,), object_sl (n_obj,), lut (n_bowtie, n_obj), slope (n_bowtie,)
    """
    spec = bh.filtered_spectrum()
    n_e = len(spec)
    kV = bh.kV[:n_e]
    bowtie_sl = np.linspace(bowtie_lengths.min(), bowtie_lengths.max(), n_bowtie)
    spec_lut = np.exp(-np.outer(bowtie_sl, bh.bowtie_ac[:n_e])) * spec  # (n_bowtie, n_e)

    object_sl = np.linspace(0, bh.object_thickness_max, int(np.ceil(bh.object_thickness_max)) * 10)
    mu = bh.object_ac[:n_e]
    # energy deposited in the scintillator, per energy bin, per unit fluence
    weight = kV * bh.scintillator_ac[:n_e] * SCINTILLATOR_THICKNESS
    fluence_0 = spec_lut @ weight  # (n_bowtie,)
    fluence_1 = spec_lut @ (np.exp(-np.outer(mu, object_sl)) * weight[:, None])  # (n_bowtie, n_obj)
    lut = -np.log(fluence_1 / fluence_0[:, None])

    # Linear model y = a*x fitted to the first 20 % of the thicknesses
    n_fit = int(round(len(object_sl) * 0.2))
    x = object_sl[:n_fit]
    slope = lut[:, :n_fit] @ x / (x @ x)
    return bowtie_sl, object_sl, lut, slope


def _inpaint_nans(img: NDArray) -> NDArray:
    """Fills NaNs smoothly: least-squares fit that makes the second differences along rows and
    columns as small as possible around them.

    Like ``inpaint_nans`` in MATLAB: smooth interpolation inside holes and linear extrapolation
    at the image edges.
    """
    nan = np.isnan(img)
    if not nan.any():
        return img
    if nan.all():
        return np.zeros_like(img)
    out = img.copy()
    idx = -np.ones(img.shape, dtype=np.int64)
    n = int(nan.sum())
    idx[nan] = np.arange(n)
    known = np.where(nan, 0.0, img)
    rows, cols, vals, rhs = [], [], [], []
    n_eq = 0
    for axis in (0, 1):
        if img.shape[axis] < 3:
            continue
        # every centred second difference (u[i-1] - 2 u[i] + u[i+1]) that touches a NaN
        sl = [slice(None)] * 2
        trip = []
        for k in range(3):
            sl[axis] = slice(k, img.shape[axis] - 2 + k)
            trip.append((idx[tuple(sl)], known[tuple(sl)]))
        touches = (trip[0][0] >= 0) | (trip[1][0] >= 0) | (trip[2][0] >= 0)
        m = int(touches.sum())
        eqs = n_eq + np.arange(m)
        b = np.zeros(m)
        for (ix, val), coef in zip(trip, (1.0, -2.0, 1.0)):
            ix, val = ix[touches], val[touches]
            unk = ix >= 0
            rows.append(eqs[unk]); cols.append(ix[unk]); vals.append(np.full(unk.sum(), coef))
            b -= np.where(unk, 0.0, coef * val)
        rhs.append(b)
        n_eq += m
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n_eq, n)).tocsc()
    out[nan] = spsolve((A.T @ A).tocsc(), A.T @ np.concatenate(rhs))
    return out


def correct_beam_hardening(log_projs: NDArray, geometry: Geometry, filepath: PathLike, scan_params) -> NDArray:
    """Beam-hardening correction of log-normalised Varian projections (BHCorrection.m).

    Args:
        log_projs: log-normalised projections, shape (N, nv, nu)
        geometry: the scan geometry
        filepath: the scan folder (with Calibrations/ASC/Factory/Calibration.xml)
        scan_params: the scan's ScanParams (for the voltage, kV filter and bowtie)

    Returns:
        NDArray: corrected projections, same shape, float32
    """
    print("Beam hardening correction (slow, and only worth it if the scan shows cupping):")
    acq_text = lambda field: scan_params._get_field(field).text  # noqa: E731
    bh = BHParams(filepath, acq_text("Voltage"), acq_text("KVFilter"), acq_text("Bowtie"))
    lengths = bowtie_path_lengths(geometry, bh)
    bowtie_sl, object_sl, lut, slope = water_lookup(bh, lengths)

    # Per pixel: the LUT for its bowtie length, inverted to an equivalent water thickness,
    # times the slope of the linear model.
    n_v, n_u = lengths.shape
    by_pixel = np.ascontiguousarray(np.moveaxis(log_projs, 0, -1))  # (nv, nu, N): one pixel's angles together
    out = np.empty(by_pixel.shape, dtype=np.float64)
    pos = np.interp(lengths, bowtie_sl, np.arange(len(bowtie_sl)))
    k0 = np.minimum(np.floor(pos).astype(int), len(bowtie_sl) - 2)
    w = pos - k0
    pix_slope = np.interp(lengths, bowtie_sl, slope)
    for iv in tqdm(range(n_v)):
        for iu in range(n_u):
            k, wk = k0[iv, iu], w[iv, iu]
            signal = (1 - wk) * lut[k] + wk * lut[k + 1]
            thickness = np.interp(by_pixel[iv, iu], signal, object_sl, left=np.nan, right=np.nan)
            out[iv, iu] = pix_slope[iv, iu] * thickness
    del by_pixel
    out = np.moveaxis(out, -1, 0)

    out[out < 0] = np.nan
    for i in range(out.shape[0]):
        out[i] = _inpaint_nans(out[i])
    out[out < 0] = np.finfo(np.float32).eps
    print("Beam hardening correction is done.")
    return out.astype(np.float32)
