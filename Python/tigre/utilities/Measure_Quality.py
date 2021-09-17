"""
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD.
%                     See the full license at
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           MATLAB (original): Manasavee Lohvithee
                      Python: Reuben Lindroos
%--------------------------------------------------------------------------
"""
from __future__ import division

from functools import reduce

import numpy as np


def Measure_Quality(res_prev, res, QualMeasOpts):  # noqa: N803
    """
    :param res_prev: (np.ndarray)
        object being compared to
    :param res: (np.ndarray)
        true image
    :param QualMeasOpts: (str or listof(str))
        any of: 'RMSE', 'nRMSE' , 'CC', 'MSSIM', 'UQI', 'SSD'
    :return:
    """
    values = []
    # make it a list:
    for qual in QualMeasOpts:
        if "RMSE" == qual:
            N = np.prod(res_prev.shape)
            diff = res_prev - res
            values.append(np.sqrt(np.sum(diff ** 2) / N))
        if "nRMSE" == qual:
            N = reduce(lambda x, y: x * y, res_prev.shape)
            diff = res_prev - res
            values.append(
                (np.sqrt(np.sum(diff ** 2) / N) / ((0.00001 + np.sqrt(np.sum(res ** 2)) / N)))
            )
        if "CC" == qual:
            cc = np.corrcoef(res_prev.ravel(), res.ravel())
            values.append(cc[0, 1])

        if "MSSIM" == qual:
            N = np.prod(res_prev.shape)

            # Compute the mean pixel values of the two images

            mean_res_p = res_prev.mean()
            mean_res = res.mean()
            if mean_res == 0 and mean_res_p == 0:
                raise ValueError("Initialising with 0 matrix not valid")
            # Luminance Comparison

            K1 = 0.01  # K1 is a small constant <<1
            d = np.max(res_prev) - np.min(res_prev)  # dynamic range of the pixel values
            l = ((2 * mean_res * mean_res_p) + (K1 * d) ** 2) / (
                (mean_res_p ** 2) + (mean_res ** 2) + K1 * d ** 2
            )

            # Contrast comparison

            K2 = 0.02
            sres_p = res_prev.std()
            sres = res.std()

            c = ((2 * sres_p * sres) + (K2 * d) ** 2) / ((sres_p ** 2) + (sres ** 2) + K2 * d ** 2)

            # Structure comparison
            diffres_p = res_prev - mean_res_p
            diffres = res - mean_res
            delta = (1 / (N - 1)) * np.sum(diffres_p * diffres)
            s = (delta + (((K2 * d) ** 2)) / 2) / ((sres_p * sres) + ((K2 * d ** 2) / 2))

            values.append((1 / N) * l * c * s)
        if "UQI" == qual:
            N = np.prod(res_prev.shape)
            # Mean
            mean_res_p = np.mean(res_prev)
            mean_res = np.mean(res)

            # Variance
            varres_p = np.var(res_prev)
            varres = np.var(res)
            if mean_res == 0 and mean_res_p == 0:
                raise ValueError("Initialising with 0 matrix not valid")
            # Covariance
            cova = np.sum((res - mean_res) * (res_prev - mean_res_p)) / (N - 1)
            front = (2 * cova) / (varres + varres_p)
            back = (2 * mean_res * mean_res_p) / ((mean_res ** 2) + (mean_res_p ** 2))

            values.append(np.sum(front * back))
        if "SSD" == qual:
            values.append(np.sum((res_prev - res) ** 2))
    if len(values) == 1:
        return np.array(values[0])
    else:
        return np.array(values)
