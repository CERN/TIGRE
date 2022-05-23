function [proj_lg, geo, angles] = VarianDataLoader(datafolder, varargin)
% VarianDataLoader: Loads Varian CBCT projection, geomtry and angles data
%   Optional parameter: Motion lag correction. Default True. 
%   gpuids: Usable GPUs. Default, the first one.
% Load all dataset that are needed for reconstruction
% Tested on TrueBeam 2.0 and 2.7
% Date: 2021-04-02
% Author: Yi Du (yi.du@hotmail.com)
% datafolder = '~/your_data_path/varian/2020-01-01_123456/';


%% Input Parser
% ACDC: acceleration & deceleration correction (default: true)
% DPS: Detector Point Spread correction (default: true)
% SC: Scatter Correction (default: true)
% BH: Beam Hardening correction (default: false, due to Bowtie, BH correction is not required at all.)
[tag_ACDC, tag_DPS, tag_SC, tag_BH, gpuids] = parse_inputs(varargin{:});

%% GPU initialization
reset(gpuDevice(gpuids.devices(0)+1));

%% Load geometry
[geo, ScanXML] = GeometryFromXML(datafolder);

%% Load Scatter Correction Calibration Parameter Structure
if(tag_DPS || tag_SC)
    ScCalib = ScCalibFromXML(datafolder);
end

%% Remove over-sampled projections due to acceleration and deceleration
thd = 0;
if(tag_ACDC)
    angular_interval = str2double(ScanXML.Acquisitions.Velocity.Text)...
        ./str2double(ScanXML.Acquisitions.FrameRate.Text);
    thd = angular_interval *0.9;
end

%% Load proj and angles
disp('Loading Proj: ')
[proj, angles, airnorm] = ProjLoader(datafolder, thd);
% Detector point scatter correction
if(tag_DPS)
    disp('Proj DPS: ')
    proj = DetectorPointScatterCorrection(proj, geo, ScCalib, gpuids);
end

%% Load blank scan
disp('Loading Blk: ')
[Blk, Sec, BlkAirNorm] = BlkLoader(datafolder);
% Detector point scatter correction
if(tag_DPS)
    disp('Blk DPS: ')
    Blk = DetectorPointScatterCorrection(Blk, geo, ScCalib);
end
%% Scatter Correction
if(tag_SC)
    disp('Scatter correction onging: ')
    proj = ScatterCorrection(ScCalib, Blk, BlkAirNorm, proj, airnorm, geo);
    disp('Scatter correction is completed.')
end

%% Airnorm and Logarithmic Normalization
proj_lg = LogNormal(proj, angles, airnorm, Blk, Sec, BlkAirNorm, gpuids);
disp('Log Normalization is completed.')
% remove anomolies
proj_lg = EnforcePositive(proj_lg); 

%% Beam Hardening correction is applied (kind of slow)
if(tag_BH)
    [proj_lg, ~] = BHCorrection(datafolder, geo, ScanXML, proj_lg, gpuids);
end

%% Remove anomalies
proj_lg = EnforcePositive(proj_lg); 

%% mediant filtering along colume-order
proj_lg = RingRemoval(proj_lg);

%% double to single
proj_lg = single(proj_lg);
angles = deg2rad(angles);

%
disp('Data processing is complete! Ready for reconstruction: ')

end


function [tag_ACDC, tag_DPS, tag_SC, tag_BH, gpuids] = parse_inputs(varargin)
% create input parser
p = inputParser;
% add optional parameters
addParameter(p,'acdc', true);
addParameter(p,'dps', true);
addParameter(p,'sc', true);
% BH performance is very lame for unclear reason
addParameter(p,'bh', false);

addParameters(p,'gpuids',GpuIds())

%execute
parse(p,varargin{:});
%extract
tag_ACDC=p.Results.acdc;
tag_DPS=p.Results.dps;
tag_SC=p.Results.sc;
tag_BH=p.Results.bh;
gpuids=p.Results.gpuids;
end
