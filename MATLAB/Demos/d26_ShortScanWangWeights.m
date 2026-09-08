%% Demo 26: Wang detector-offset weights on a short scan
%
% Wang weighting compensates a displaced detector on a FULL circular scan: it
% ramps one side of the detector down and relies on the opposing (beta + pi)
% views to bring the coverage back to uniform. On a short scan those opposing
% views do not exist, so the ramp survives into the reconstruction as a
% one-sided shading (a "teardrop"). Any non-zero offDetector switches the
% weights on - a sub-pixel calibrated offset is enough - so a short scan with a
% calibrated geometry was silently shaded.
%
% FDK now skips the Wang weights when the angles do not cover a full circle
% (the same test it already uses to decide Parker weighting). This demo
% reconstructs a uniform cylinder from a 200-degree scan three ways and prints
% the left/right interior ratio:
%
%   centred detector                    -> balanced (reference)
%   0.1 mm offset, Wang applied by hand -> the teardrop (what happened before)
%   0.1 mm offset, FDK default          -> Wang skipped, balanced again
%
% Wang weights are still applied on full scans, where they belong.
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD.
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
%--------------------------------------------------------------------------
clear; close all;

%% Geometry: default, 64^3 image
geo = defaultGeometry('nVoxel', [64; 64; 64]);
geo.sVoxel = [64; 64; 64];
geo.dVoxel = geo.sVoxel ./ geo.nVoxel;

%% Uniform cylinder, mu = 0.02
[xx, yy, zz] = meshgrid(1:64, 1:64, 1:64);
mu = 0.02;
phantom = single(((xx - 32.5).^2 + (yy - 32.5).^2 <= 20^2) * mu);
core = (xx - 32.5).^2 + (yy - 32.5).^2 <= 14^2;

%% Short scan: 200 degrees
angles = linspace(0, deg2rad(200), 200);
proj = Ax(phantom, geo, angles);

lr = @(vol) [mean(vol(core & xx <= 32 & zz >= 17 & zz <= 48)), ...
             mean(vol(core & xx >  32 & zz >= 17 & zz <= 48))];

%% 1. centred detector (reference)
vol_ref = FDK(proj, geo, angles);
r = lr(vol_ref);
fprintf('centred detector           : left %.5f right %.5f ratio %.3f\n', r(1), r(2), r(1) / r(2));

%% 2. 0.1 mm offset, Wang weights applied regardless of the arc (the old behaviour)
geo_off = geo;
geo_off.offDetector = [0.1; 0];          % u offset, a fraction of a pixel
proj_w = proj .* redundancy_weighting(geo_off);
vol_wang = FDK(proj_w, geo_off, angles, 'wang', false);
r = lr(vol_wang);
fprintf('offset, Wang applied (old) : left %.5f right %.5f ratio %.3f\n', r(1), r(2), r(1) / r(2));

%% 3. 0.1 mm offset, FDK default: Wang skipped on the short scan (with a warning)
vol_gated = FDK(proj, geo_off, angles);
r = lr(vol_gated);
fprintf('offset, FDK default (new)  : left %.5f right %.5f ratio %.3f\n', r(1), r(2), r(1) / r(2));

%% Show the central axial slice of each
figure('Name', '200-degree scan: Wang detector-offset weights need a full circle');
subplot(1, 3, 1); imshow(squeeze(vol_ref(:, :, 32)),   [0 1.5 * mu]); title('centred (reference)');
subplot(1, 3, 2); imshow(squeeze(vol_wang(:, :, 32)),  [0 1.5 * mu]); title('offset, Wang applied (before)');
subplot(1, 3, 3); imshow(squeeze(vol_gated(:, :, 32)), [0 1.5 * mu]); title('offset, Wang skipped (after)');
