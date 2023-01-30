function proj = DetectorPointScatterCorrection(proj, geo, ScCalib, gpuids)
%% Detector Point Scatter Correction
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Date: 2021-03-26
% Author: Yi Du (yi.du@hotmail.com)


%% Empirical values from reference paper
% unit: cm-2
% a0 = 3.43; (in refrence paper, but incorrect)
a0 = 1;

a1 = 0.000309703536035;
a1 = str2double(ScCalib.CalibrationResults.Globals.DetectorScatterModel.PScFit0.Text);

% unit: cm-1
a2 = 0.546566915157327;
a2 = str2double(ScCalib.CalibrationResults.Globals.DetectorScatterModel.PScFit1.Text);

% unit: cm-1
a3 = 0.311272841141691;
a3 = str2double(ScCalib.CalibrationResults.Globals.DetectorScatterModel.PScFit2.Text);

% unit: cm-1
a4 = 0.002472148007134;
a4 = str2double(ScCalib.CalibrationResults.Globals.DetectorScatterModel.PScFit3.Text);

a5 = -12.6606856375944;
a5 = str2double(ScCalib.CalibrationResults.Globals.DetectorScatterModel.PScFit4.Text);

% for amplitude normalization
CoverSPR = 0.04;

% unit: mm
offset=geo.offDetector;

% grid unit: mm
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1); % + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2); % + offset(2);
% unit mm - > cm
% unit converter: 1/10 for mm-> cm
mm2cm = 1/10;
us = us * mm2cm;
vs = vs * mm2cm;

%% Downsampling
% about 10 mm in axial direction
dus = downsample(us, 26);
% about 4 mm in axial direction
dvs = downsample(vs, 10);

ds_rate = 8;
dus = decimate(us, ds_rate);
dvs = decimate(vs, ds_rate);


%% Grid mesh
[uu,vv] = meshgrid(us,vs); %detector
[duu, dvv] = meshgrid(dus,dvs); %detector

%% Scatter convolution kernel
grid = sqrt(duu.^2 + dvv.^2);
% a0 is the normalization factor
hd = a0*( a1* exp(-a2 * grid) + a3 * (exp( -a4 * ( grid - a5).^3 )));

% a0 = CoverSPR/sum(hd(:));

% normalized to 4% SPR coverage
hd = CoverSPR/sum(hd(:)) .* hd;

%% GPU based
% reset(gpuDevice(1));
reset(gpuDevice(gpuids.devices(0)+1));
gproj = gpuArray(single(proj));

%% 2D Convolution with downsampling and upsampling
for ii = 1:size(proj,3)
    % CPU version
    %{
    page = interp2(uu, vv, proj(:,:,ii), duu, dvv);
    % gross scatter distribution
    sc = conv2(page, hd, 'same');
    % upsample the scatter distribution to the same grid level as the
    % measured intensity
    scpage = interp2(duu, dvv, sc, uu, vv, 'spline');
    % primary = measure - scatter
    proj(:,:,ii) = proj(:,:,ii) - scpage;
    %}
    
    % GPU version
    page = interp2(uu, vv, gproj(:,:,ii), duu, dvv);
    % gross scatter distribution
    sc = gather(conv2(page, hd, 'same'));
    % upsample the scatter distribution to the same grid level as the
    % measured intensity
    scpage = interp2(duu, dvv, sc, uu, vv, 'spline');
    % primary = measure - scatter
    gproj(:,:,ii) = gproj(:,:,ii) - scpage;
end

proj = single(gather(gproj));
% Reset GPU
reset(gpuDevice(gpuids.devices(0)+1));


%% Cutoff for over-correction
proj(proj<0) = NaN;
for ii = 1:size(proj,3)
    % default fast extrapolation: robust for noise and holes at boundaries
    proj(:,:,ii) = single(inpaint_nans(double(proj(:,:,ii)), 2));
end
proj(proj<0) = eps;

end
