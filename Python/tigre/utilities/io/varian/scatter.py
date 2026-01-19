import numpy as np
import glob
import os
import xml.etree.ElementTree as ET
from scipy.ndimage import gaussian_filter, median_filter
from scipy.fft import fft2, ifft2
from scipy.signal import decimate, convolve2d
from scipy.interpolate import interpn
from Python.tigre.utilities.io.varian.utils import mm2cm, XML
from tqdm import tqdm


class DetScattParams(XML):
    def __init__(self, filepath, xml_reader):
        super().__init__(filepath, xml_reader)
        self.params = self._get_params()

    def _get_params(self):
        det_scatter_model = self.root.find(
            "CalibrationResults/Globals/DetectorScatterModel", self.ns
        )
        return [float(elem.text) for elem in det_scatter_model]

    def calculate_dps_kernel(self, grid_coords):
        U, V = grid_coords
        grid = np.sqrt(U**2 + V**2)
        det_kernel = self.params[0] * np.exp(-self.params[1] * grid) + self.params[2] * (
            np.exp(-self.params[3] * (grid - self.params[4]) ** 3)
        )
        det_kernel = self.params[5] * det_kernel / np.sum(det_kernel)
        return det_kernel


class ScattParams(XML):
    def __init__(self, filepath, xml_reader):
        super().__init__(filepath, xml_reader)
        self.obj_scatt_models = self._get_obj_scatt_models()
        self.obj_scatt_fits = self._get_obj_scatt_fits()
        self.thickness_params = self._get_thickness_params()
        self.grid_efficiency = self._get_grid_efficiency()
        self.gamma = self._get_prop_coeff()
        self.form_func_params = self._get_form_func_params()
        self.ampl_params = self._get_ampl_params()

    def _get_obj_scatt_models(self):
        return self.root.find("CalibrationResults/ObjectScatterModels", self.ns)

    def _get_obj_scatt_fits(self):
        return self.obj_scatt_models.findall("ObjectScatterModel/ObjectScatterFit", self.ns)

    def _get_thickness_params(self):
        mu_water = float(self.root.find("CalibrationResults/Globals/muH2O", self.ns).text)  # mm^-1
        sigma_u = float(self.root.find("CalibrationResults/Globals/AsymPertSigmaMMu", self.ns).text)
        sigma_v = float(self.root.find("CalibrationResults/Globals/AsymPertSigmaMMv", self.ns).text)
        bounds = [
            float(thickness.text)
            for thickness in self.obj_scatt_models.findall("ObjectScatterModel/Thickness", self.ns)
        ]
        thickness_params = dict(mu_water=mu_water, sigma_u=sigma_u, sigma_v=sigma_v, bounds=bounds)
        return thickness_params

    def _get_grid_efficiency(self):
        grid_efficiency = float(
            self.obj_scatt_models.find(
                "ObjectScatterModel/GridEfficiency/LamellaTransmission", self.ns
            ).text
        )
        return grid_efficiency

    def _get_prop_coeff(self):
        gamma = float(
            self.obj_scatt_models.find("ObjectScatterModel/ObjectScatterFit/gamma", self.ns).text
        )
        return gamma

    def _get_form_func_params(self):
        sigma1 = [
            float(elem.find("sigma1", self.ns).text) for elem in self.obj_scatt_fits
        ]  # (cm^-1)
        sigma2 = [
            float(elem.find("sigma2", self.ns).text) for elem in self.obj_scatt_fits
        ]  # (cm^-1)
        B = [float(elem.find("B", self.ns).text) for elem in self.obj_scatt_fits]
        return dict(sigma1=sigma1, sigma2=sigma2, B=B)

    def _get_ampl_params(self):
        A = [float(elem.find("A", self.ns).text) for elem in self.obj_scatt_fits]
        alpha = [float(elem.find("alpha", self.ns).text) for elem in self.obj_scatt_fits]
        beta = [float(elem.find("beta", self.ns).text) for elem in self.obj_scatt_fits]
        return dict(A=A, alpha=alpha, beta=beta)

    def num_thicknesses(self):
        return len(self.thickness_params["bounds"])

    def create_thickness_masks(self, thickness_map):
        mask = np.zeros_like(thickness_map)
        masks = np.tile(mask[np.newaxis, :, :], [self.num_thicknesses(), 1, 1])
        for i in range(self.num_thicknesses() - 1):
            masks[i] = (thickness_map > self.thickness_params["bounds"][i]) * (
                thickness_map < self.thickness_params["bounds"][i + 1]
            )
        masks[-1] = thickness_map > self.thickness_params["bounds"][-1]
        return masks

    def calculate_grid_response(self, u, v):
        # TODO: read params from Calibration.xml
        K = -0.15
        B = 1
        t_ratio = K * abs(v.T) + B
        kernel = np.tile(t_ratio[:, np.newaxis], [1, len(u)])
        kernel[kernel < self.grid_efficiency] = self.grid_efficiency
        return kernel

    def estimate_water_equiv_thickness(self, blank, proj, smooth=False, sigma=None):
        log_norm = _log_normalize(blank, proj)
        thickness = log_norm / self.thickness_params["mu_water"]
        if smooth:
            if sigma is None:
                raise ValueError("sigma cannot be 'None' if smooth=True.")
            thickness = gaussian_filter(thickness, sigma)
        return thickness

    def calculate_form_functions(self, thickness_map, grid_coords):
        U, V = grid_coords
        grid = U**2 + V**2
        gform = np.zeros_like(thickness_map)
        gforms = np.tile(gform[np.newaxis, :, :], [self.num_thicknesses(), 1, 1])

        for i in range(self.num_thicknesses()):
            gforms[i] = np.exp(
                -0.5 * grid / (self.form_func_params["sigma1"][i] ** 2)
            ) + self.form_func_params["B"][i] * np.exp(
                -0.5 * grid / (self.form_func_params["sigma2"][i] ** 2)
            )
        return gforms

    def calculate_amplitudes(self, blank, proj, edge_weight):
        log_norm = _log_normalize(blank, proj)
        eps = np.finfo(proj.dtype).eps
        norm = (proj + eps) / (blank + eps)
        norm[norm > 1] = 1
        norm[norm < 0] = 0
        amplitudes = []
        for i in range(self.num_thicknesses()):
            amplitude = (
                mm2cm(self.ampl_params["A"][i])
                * edge_weight
                * (norm ** self.ampl_params["alpha"][i])
                * (log_norm ** self.ampl_params["beta"][i])
            )
            amplitudes.append(amplitude)
        return np.array(amplitudes)

    def calculate_thickness_sigma(self, step_du, step_dv):
        sigma = (
            mm2cm(self.thickness_params["sigma_v"]) / step_dv,
            mm2cm(self.thickness_params["sigma_u"]) / step_du,
        )
        return sigma

    @staticmethod
    def calculate_edge_response(thickness_map, threshold=50, filter_size=25, num_iter=5):
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


def read_scatt_calib_xml(filepath):
    file = glob.glob(os.path.join(filepath, "Calibrations", "SC-*", "Factory", "Calibration.xml"))
    if len(file) == 0:
        raise RuntimeError("Scatter calibration file not found.")
    elif len(file) > 1:
        raise RuntimeError("Multiple calibration files found.")
    else:
        tree = ET.parse(file[0])
        sc_calib = tree.getroot()
        return sc_calib


def _log_normalize(blank, proj):
    eps = np.finfo(proj.dtype).eps
    ratio = (blank + eps) / (proj + eps)
    ratio[ratio < 1] = 1
    log_proj = np.log(ratio)
    return log_proj


def _rescale(x, min_val=0.0, max_val=1.0):
    x_rescale = (max_val - min_val) * (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x)) + min_val
    return x_rescale


def _update_scatter_estimate(
    proj, thickness_map, thickness_masks, amplitudes, gforms, grid_kernel, gamma
):
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


def _update_primary_estimate(primary, scatter_old, scatter, lam=0.6):
    primary += lam * (scatter_old - scatter)
    eps = np.finfo(primary.dtype).eps
    primary[primary < eps] = eps
    return primary


def _calculate_primary(proj, scatt, max_scatt_frac=0.95):
    scatt_frac = scatt / (proj + np.finfo(proj.dtype).eps)

    scatt_frac = median_filter(scatt_frac, size=3)
    scatt_frac[scatt_frac < 0] = 0
    scatt_frac = np.minimum(scatt_frac, max_scatt_frac)
    return proj * (1 - scatt_frac)


def _get_detector_coords(geometry, downsample=0):
    # Detector coords,centered (cm)
    u = (
        0.5
        * (geometry.sDetector[1] - geometry.dDetector[1])
        * np.linspace(-1, 1, int(geometry.nDetector[1]))
    )
    v = (
        0.5
        * (geometry.sDetector[0] - geometry.dDetector[0])
        * np.linspace(-1, 1, int(geometry.nDetector[0]))
    )

    if downsample:
        u = decimate(u, downsample)
        v = decimate(v, downsample)

    return mm2cm(u), mm2cm(v)


def correct_detector_scatter(projs, geometry, dps_calib, downsample=8):
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
    proj_data,
    blank_proj_data,
    geometry,
    sc_calib,
    downsample=12,
    num_iter=8,
    lam=0.005,
):

    u, v = _get_detector_coords(geometry)
    U, V = np.meshgrid(u, v)

    du, dv = _get_detector_coords(geometry, downsample=downsample)
    DU, DV = np.meshgrid(du, dv)

    step_du = np.mean(np.diff(du))
    step_dv = np.mean(np.diff(dv))
    sigma = sc_calib.calculate_thickness_sigma(step_du, step_dv)

    grid_kernel = sc_calib.calculate_grid_response(du, dv)

    primaries = np.zeros_like(proj_data.projs)
    print("Performing ASKS scatter correction: ")
    for i, proj in tqdm(enumerate(proj_data.projs)):
        blank, blank_airnorm = blank_proj_data.interp_proj(proj_data.angles[i])

        cf_air = proj_data.airnorms[i] / blank_airnorm
        blank = interpn((v, u), blank * cf_air, (DV, DU))
        blank = np.float32(blank)
        primary = interpn((v, u), proj, (DV, DU))
        primary = np.float32(primary)
        scatter = np.zeros_like(primary)

        for n in range(num_iter):
            scatter_old = scatter
            thickness_map = sc_calib.estimate_water_equiv_thickness(
                blank, primary, smooth=True, sigma=sigma
            )
            edge_weight = sc_calib.calculate_edge_response(thickness_map)
            thickness_masks = sc_calib.create_thickness_masks(thickness_map)
            gforms = sc_calib.calculate_form_functions(thickness_map, (DU, DV))
            amplitudes = sc_calib.calculate_amplitudes(blank, primary, edge_weight)

            scatter = _update_scatter_estimate(
                primary,
                mm2cm(thickness_map),
                thickness_masks,
                amplitudes,
                gforms,
                grid_kernel,
                sc_calib.gamma,
            )
            primary = _update_primary_estimate(primary, scatter_old, scatter, lam=lam)

        # Upsample
        scatter_est = interpn(
            (dv, du), scatter, (V, U), method="cubic", bounds_error=False, fill_value=None
        )
        primaries[i] = _calculate_primary(proj, scatter_est)
    return primaries
