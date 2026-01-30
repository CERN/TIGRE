from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
import glob
import os
import xml.etree.ElementTree as ET
from scipy.ndimage import gaussian_filter, median_filter
from scipy.fft import fft2, ifft2
from scipy.signal import decimate, convolve2d
from scipy.interpolate import interpn
from tigre.utilities.geometry import Geometry
from tigre.utilities.io.varian.utils import cm2mm, XML, PathLike, XMLReader
from tigre.utilities.io.varian.varian_io import ProjData
from tqdm import tqdm


def _read_scatt_xml(filepath: PathLike) -> ET.Element:
    """Reads scatter correction parameters from Calibration.xml stored in
      '/Calibrations/SC-*/Factory' directory. Raises runtime error if none or multiple files found.

    Args:
        filepath (PathLike): scan folder

    Returns:
        ET.Element: xml root element
    """
    file = glob.glob(os.path.join(filepath, "Calibrations", "SC-*", "Factory", "Calibration.xml"))
    if len(file) == 0:
        raise RuntimeError("Scatter calibration file not found.")
    elif len(file) > 1:
        raise RuntimeError("Multiple calibration files found.")
    else:
        tree = ET.parse(file[0])
        sc_calib = tree.getroot()
        return sc_calib


class DetScattParams(XML):
    """Class for detector scatter correction parameters (read from Calibration.xml file)."""

    def __init__(self, filepath: PathLike, xml_reader: XMLReader = _read_scatt_xml) -> None:
        super().__init__(filepath, xml_reader)
        self.params: list[float] = self._get_params()

    def _get_params(self) -> list[float]:
        dps_elems = self._get_field("CalibrationResults/Globals/DetectorScatterModel")
        try:
            dps_params = [float(elem.text) for elem in dps_elems]
        except Exception:
            raise ValueError("Invalid DetectorScatterModel parameters.")
        else:
            dps_params[1] /= cm2mm(1)  # convert to mm^-1
            dps_params[3] /= cm2mm(1)  # convert to mm^-1
            return dps_params

    def calculate_dps_kernel(self, grid_coords: tuple[NDArray, NDArray]) -> NDArray:
        U, V = grid_coords
        grid = np.sqrt(U**2 + V**2)

        det_kernel = self.params[0] * np.exp(-self.params[1] * grid) + self.params[2] * (
            np.exp(-self.params[3] * (grid - self.params[4]) ** 3)
        )
        det_kernel = self.params[5] * det_kernel / np.sum(det_kernel)
        return det_kernel


class ScattParams(XML):
    """Class for FASKS scatter correction parameters (read from Calibration.xml file)."""

    def __init__(self, filepath: PathLike, xml_reader: XMLReader = _read_scatt_xml) -> None:
        super().__init__(filepath, xml_reader)
        self.obj_scatt_models: ET.Element = self._get_obj_scatt_models()
        self.obj_scatt_fits: list[ET.Element] = self._get_obj_scatt_fits()
        self.thickness_params: dict[str, float | list[float]] = self._get_thickness_params()
        self.grid_efficiency: float = self._get_grid_efficiency()
        self.gamma: float = self._get_prop_coeff()
        self.form_func_params: dict[str, NDArray] = self._get_form_func_params()
        self.ampl_params: dict[str, NDArray] = self._get_ampl_params()

    def _get_obj_scatt_models(self) -> ET.Element:
        return self._get_field(
            "CalibrationResults/ObjectScatterModels",
        )

    def _get_obj_scatt_fits(self) -> list[ET.Element]:
        return self.obj_scatt_models.findall("ObjectScatterModel/ObjectScatterFit", self.ns)

    def _get_thickness_params(self) -> dict[str, float | list[float]]:
        mu_water = float(self._get_field("CalibrationResults/Globals/muH2O").text)  # mm^-1
        sigma_u = float(self._get_field("CalibrationResults/Globals/AsymPertSigmaMMu").text)  # mm
        sigma_v = float(self._get_field("CalibrationResults/Globals/AsymPertSigmaMMv").text)  # mm
        bounds = [
            float(thickness.text)
            for thickness in self.obj_scatt_models.findall("ObjectScatterModel/Thickness", self.ns)
        ]  # mm
        thickness_params = dict(mu_water=mu_water, sigma_u=sigma_u, sigma_v=sigma_v, bounds=bounds)
        return thickness_params

    def _get_grid_efficiency(self) -> float:
        grid_efficiency = float(
            self.obj_scatt_models.find(
                "ObjectScatterModel/GridEfficiency/LamellaTransmission", self.ns
            ).text
        )
        return grid_efficiency

    def _get_prop_coeff(self) -> float:
        gamma = float(
            self.obj_scatt_models.find("ObjectScatterModel/ObjectScatterFit/gamma", self.ns).text
        )
        return gamma

    def _get_form_func_params(self) -> dict[str, NDArray]:
        sigma1 = np.array(
            [float(elem.find("sigma1", self.ns).text) for elem in self.obj_scatt_fits]
        )  # cm
        sigma2 = np.array(
            [float(elem.find("sigma2", self.ns).text) for elem in self.obj_scatt_fits]
        )  # cm
        B = np.array([float(elem.find("B", self.ns).text) for elem in self.obj_scatt_fits])

        sigma1 = cm2mm(sigma1)  # convert to mm
        sigma2 = cm2mm(sigma2)  # convert to mm

        return dict(sigma1=sigma1, sigma2=sigma2, B=B)

    def _get_ampl_params(self) -> dict[str, NDArray]:
        A = np.array([float(elem.find("A", self.ns).text) for elem in self.obj_scatt_fits])  # cm^-2
        alpha = np.array([float(elem.find("alpha", self.ns).text) for elem in self.obj_scatt_fits])
        beta = np.array([float(elem.find("beta", self.ns).text) for elem in self.obj_scatt_fits])
        A /= cm2mm(1) ** 2  # convert to mm^-2
        return dict(A=A, alpha=alpha, beta=beta)

    def num_thicknesses(self) -> int:
        return len(self.thickness_params["bounds"])

    def create_thickness_masks(self, thickness_map: NDArray) -> NDArray:
        mask = np.zeros_like(thickness_map)
        masks = np.tile(mask[np.newaxis, :, :], [self.num_thicknesses(), 1, 1])
        for i in range(self.num_thicknesses() - 1):
            masks[i] = (thickness_map > self.thickness_params["bounds"][i]) * (
                thickness_map < self.thickness_params["bounds"][i + 1]
            )
        masks[-1] = thickness_map > self.thickness_params["bounds"][-1]
        return masks

    def calculate_grid_response(self, u: NDArray, v: NDArray) -> NDArray:
        # TODO: read params from Calibration.xml
        K = -0.15 / cm2mm(1)
        B = 1
        t_ratio = K * np.abs(v.T) + B
        kernel = np.tile(t_ratio[:, np.newaxis], [1, len(u)])
        kernel[kernel < self.grid_efficiency] = self.grid_efficiency
        return kernel

    def estimate_water_equiv_thickness(
        self,
        blank: NDArray,
        proj: NDArray,
        pixel_size: tuple[float, float],
    ) -> NDArray:
        log_norm = _log_normalize(blank, proj)
        thickness_map = log_norm / self.thickness_params["mu_water"]
        sigma = self.calculate_thickness_sigma(*pixel_size)
        thickness_map = gaussian_filter(thickness_map, sigma)
        return thickness_map

    def calculate_form_functions(
        self, thickness_map: NDArray, grid_coords: tuple[NDArray, NDArray]
    ) -> NDArray:
        U, V = grid_coords
        grid = U**2 + V**2
        gform = np.zeros_like(thickness_map)
        gforms = np.tile(gform[np.newaxis, :, :], [self.num_thicknesses(), 1, 1])

        for i in range(self.num_thicknesses()):
            gforms[i] = np.exp(
                -0.5 * grid / ((self.form_func_params["sigma1"][i]) ** 2)
            ) + self.form_func_params["B"][i] * np.exp(
                -0.5 * grid / ((self.form_func_params["sigma2"][i]) ** 2)
            )
        return gforms

    def calculate_amplitudes(self, blank: NDArray, proj: NDArray, edge_weight: NDArray) -> NDArray:
        log_norm = _log_normalize(blank, proj)
        eps = np.finfo(proj.dtype).eps
        norm = (proj + eps) / (blank + eps)
        norm[norm > 1] = 1
        norm[norm < 0] = 0
        amplitudes = []
        for i in range(self.num_thicknesses()):
            amplitude = (
                (self.ampl_params["A"][i])
                * edge_weight
                * (norm ** self.ampl_params["alpha"][i])
                * (log_norm ** self.ampl_params["beta"][i])
            )
            amplitudes.append(amplitude)
        return np.array(amplitudes)

    def calculate_thickness_sigma(self, step_du: float, step_dv: float) -> tuple[float, float]:
        sigma = (
            self.thickness_params["sigma_v"] / step_dv,
            self.thickness_params["sigma_u"] / step_du,
        )
        return sigma

    @staticmethod
    def calculate_edge_response(
        thickness_map: NDArray,
        threshold: float = 50.0,
        filter_size: int = 25,
        num_iter: int = 5,
    ) -> NDArray:
        # TODO: read params from calibration.xml file
        edge_weight = thickness_map > threshold
        edge_weight = np.array(edge_weight, dtype="float")
        tmp_mask = edge_weight
        avg_kernel = np.ones([filter_size, filter_size]) / filter_size

        for i in range(num_iter):
            edge_weight = convolve2d(edge_weight, avg_kernel, mode="same")

        tmp = tmp_mask * edge_weight
        tmp_masked = np.ma.masked_not_equal(tmp, 0)
        tmp[tmp_masked.mask] = _rescale(tmp[tmp_masked.mask], min_val=0.6)
        edge_weight = tmp
        return edge_weight


def _log_normalize(blank: NDArray, proj: NDArray) -> NDArray:
    eps = np.finfo(proj.dtype).eps
    ratio = (blank + eps) / (proj + eps)
    ratio[ratio < 1] = 1
    log_proj = np.log(ratio)
    return log_proj


def _rescale(x: NDArray, min_val: float = 0.0, max_val: float = 1.0) -> NDArray:
    """Rescales input x to range [min_val, max_val] (default: [0,1])."""
    x_rescale = (max_val - min_val) * (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x)) + min_val
    return x_rescale


def _update_scatter_estimate(
    proj: NDArray,
    thickness_map: NDArray,
    thickness_masks: NDArray,
    amplitudes: NDArray,
    gforms: NDArray,
    grid_kernel: NDArray,
    gamma: float,
) -> NDArray:
    num_groups = len(thickness_masks)
    rep_proj = np.tile(proj[np.newaxis, :, :], [num_groups, 1, 1])
    term1 = rep_proj * thickness_masks * amplitudes
    rep_grid_kernel = np.tile(grid_kernel[np.newaxis, :, :], [num_groups, 1, 1])
    term2 = fft2(gforms * rep_grid_kernel)
    rep_thickness = np.tile(thickness_map[np.newaxis, :, :], [num_groups, 1, 1])
    term3 = fft2(rep_thickness * term1)
    tmp1 = np.sum(fft2(term1) * term2, axis=0)
    tmp2 = np.sum(term2 * term3, axis=0)
    itmp1 = ifft2(tmp1)
    itmp2 = ifft2(tmp2)
    scatter = (1 - gamma * thickness_map) * np.real(itmp1) + gamma * np.real(itmp2)
    return scatter


def _update_primary_estimate(
    primary: NDArray, scatter_old: NDArray, scatter: NDArray, lam: float = 0.6
) -> NDArray:
    primary += lam * (scatter_old - scatter)
    eps = np.finfo(primary.dtype).eps
    primary[primary < eps] = eps
    return primary


def _calculate_primary(proj: NDArray, scatt: NDArray, max_scatt_frac: float = 0.95) -> NDArray:
    scatt_frac = scatt / (proj + np.finfo(proj.dtype).eps)

    scatt_frac = median_filter(scatt_frac, size=3)
    scatt_frac = np.clip(scatt_frac, 0, max_scatt_frac)
    return proj * (1 - scatt_frac)


def _get_detector_coords(geometry: Geometry, downsample: int = 0) -> tuple[NDArray, NDArray]:
    """Calculates the local detector coords (mm), centered at the origin.

    Args:
        geometry (Geometry): scan geometry
        downsample (int, optional): downsample factor. Defaults to 0.

    Returns:
        tuple[NDArray, NDArray]: detector coords (u,v)->(j,i)
    """
    u = 0.5 * np.linspace(-1, 1, int(geometry.nDetector[1]))
    u *= geometry.sDetector[1] - geometry.dDetector[1]
    v = 0.5 * np.linspace(-1, 1, int(geometry.nDetector[0]))
    v *= geometry.sDetector[0] - geometry.dDetector[0]

    if downsample:
        u = decimate(u, downsample)
        v = decimate(v, downsample)

    return u, v


def correct_detector_scatter(
    projs: NDArray, geometry: Geometry, dps_calib: DetScattParams, downsample: int = 8
) -> NDArray:
    """Performs detector point scatter (dps) correction using calibration factors in dps_calib.

    Args:
        projs (NDArray): uncorrected projection data
        geometry (Geometry): scan geometry
        dps_calib (DetScattParams): Calibration parameters (from xml file)
        downsample (int, optional): downsample factor. (Default: 8).

    Returns:
        NDArray: dps-corrected projections
    """
    u, v = _get_detector_coords(geometry)
    U, V = np.meshgrid(u, v)
    du, dv = _get_detector_coords(geometry, downsample=downsample)
    DU, DV = np.meshgrid(du, dv)

    det_kernel = dps_calib.calculate_dps_kernel(grid_coords=(DU, DV))

    corrected_projs = np.zeros_like(projs)
    print("Performing detector point scatter correction: ")
    for i, proj in tqdm(enumerate(projs)):
        proj_down = interpn((v, u), projs[0], (DV, DU))
        scatter_down = convolve2d(proj_down, det_kernel, mode="same")
        scatter = interpn(
            (dv, du), scatter_down, (V, U), method="cubic", bounds_error=False, fill_value=None
        )
        corrected_projs[i] = proj - scatter
    eps = np.finfo(proj.dtype).eps
    corrected_projs[corrected_projs < eps] = eps

    return corrected_projs


def correct_scatter(
    proj_data: ProjData,
    blank_proj_data: ProjData,
    geometry: Geometry,
    sc_calib: ScattParams,
    downsample: int = 12,
    max_iter: int = 8,
    lam: float = 0.005,
    min_delta: float = 1e-16,
) -> NDArray:
    """Performs FASKS scatter correction (sc) based on algorithm described in Sun & Star-Lack 2010
    (doi: 10.1088/0031-9155/55/22/007).

    Args:
        proj_data (ProjData): uncorrected projection data (assumes dps correction already applied)
        blank_proj_data (ProjData): blank projection data (assumes dps correction already applied)
        geometry (Geometry): scan geometry
        sc_calib (ScattParams): scatter calibration parameters (from xml file)
        downsample (int, optional): downsample factor. (Default: 12).
        max_iter (int, optional): maximum number of iterations for scatter estimation. (Default: 8).
        lam (float, optional): relaxation factor. (Default: 0.005).
        min_delta (float, optional): minimum change in scatter estimation (mean abs. diff.). (Default: 1e-16).

    Returns:
        NDArray: scatter-corrected projections
    """

    u, v = _get_detector_coords(geometry)
    U, V = np.meshgrid(u, v)

    du, dv = _get_detector_coords(geometry, downsample=downsample)
    DU, DV = np.meshgrid(du, dv)

    step_du = np.mean(np.diff(du))
    step_dv = np.mean(np.diff(dv))

    grid_kernel = sc_calib.calculate_grid_response(du, dv)

    primaries = np.zeros_like(proj_data.projs)

    print("Performing FASKS scatter correction: ")
    for i, proj in tqdm(enumerate(proj_data.projs)):
        blank = blank_proj_data.interp_proj(proj_data.angles[i])
        blank = interpn((v, u), blank, (DV, DU))
        primary = interpn((v, u), proj, (DV, DU))
        scatter = np.zeros_like(primary)

        n_iter = 0
        delta_scatter = 1.0
        while delta_scatter > min_delta and n_iter < max_iter:
            scatter_old = scatter
            thickness_map = sc_calib.estimate_water_equiv_thickness(
                blank, primary, pixel_size=(step_du, step_dv)
            )
            edge_weight = sc_calib.calculate_edge_response(thickness_map)
            thickness_masks = sc_calib.create_thickness_masks(thickness_map)
            gforms = sc_calib.calculate_form_functions(thickness_map, (DU, DV))
            amplitudes = sc_calib.calculate_amplitudes(blank, primary, edge_weight)

            scatter = _update_scatter_estimate(
                primary,
                thickness_map,
                thickness_masks,
                amplitudes,
                gforms,
                grid_kernel,
                sc_calib.gamma,
            )
            primary = _update_primary_estimate(primary, scatter_old, scatter, lam=lam)
            delta_scatter = np.mean(np.abs(scatter - scatter_old))
            n_iter += 1

        scatter_est = interpn(
            (dv, du), scatter, (V, U), method="cubic", bounds_error=False, fill_value=None
        )
        primaries[i] = _calculate_primary(proj, scatter_est)
    return primaries
