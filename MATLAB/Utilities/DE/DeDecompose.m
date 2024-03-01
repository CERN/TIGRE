function [material1_proj,material2_proj, angles_out] = DeDecompose(proj_l, angles_l , proj_h, angles_h, varargin)
%DEDECOMPSE Maps CBCT projections taken at two kVp settings to equivalent
%thicknesses of material1 and material2 using precalibrated transfer functions.
%   Coefficients C and D are the coefficients for the third-degree transfer
%   functions (may need to be modified for different systems)
%
%   Other basis materials and kV settings can be employed by adding
%   additional calibration data to "de_cal.mat"
%
%   Proj_l and angles_l are the projection and angle data for the low-kV
%   scan
%
%   Proj_h and angles_h are the projection and angle data for the high-kV
%   scan
%
%   Optional Arguments
%       scanner - default VarianTrueVeam
%
%       kVl, kVh - the kV settings for the low- and high-kV scans,
%       respectively.  defaults: 80, 140
%
%       interp_model - model for sinogram interpolation
%
%           'nearest' - nearest neighbor (default).
%
%           'skip' - same as 'nearest, but skips projections without close
%           matches.  Requires setting 'tolerance' value.
%
%           'linear' - linear interpolation (unreliable for large
%           projection separations)
%
%   Coded by: Andrew Keeler
%   Contact: akeeler@luc.edu
%   Date: February 28, 2024

% Load calibrated coefficients for transfer functions


[scanner, calibration_pair, model, tolerance, material1, material2] = parse_inputs(varargin);

% read NIST database (values in TIGRE)
fid = fopen('./../../../Common/data/dual_energy_calibrations.json');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
calibration_data = jsondecode(str);

if ~isfield(calibration_data, scanner)
    error(['Scanner: ', scanner, ' has no calibration data'])
end
if ~isfield(calibration_data.(scanner),calibration_pair)
    error(['Calibration pair of energies: ', calibration_pair, ' not found in the scanner data'])
end
materials = fieldnames(calibration_data.(scanner).(calibration_pair));

if ~isfield(calibration_data.(scanner).(calibration_pair),material1)
    disp([material1, ' not a (yet) supported material in the scanner at the given energies'])
    disp("Supported materials:")
    for i=1:len(materials)
        disp(materials{i})
    end
    error("")
end
if ~isfield(calibration_data.(scanner).(calibration_pair),material2)
    disp([material2, ' not a (yet) supported material in the scanner at the given energies'])
    disp("Supported materials:")
    for i=1:len(materials)
        disp(materials{i})
    end
    error("")
end

c = calibration_data.(scanner).(calibration_pair).(material1);
d = calibration_data.(scanner).(calibration_pair).(material2);

%% check that the projections are the same size
for n = 1:2
    if size(proj_h,n) ~= size(proj_l,n)
        fprintf('High- and low-energy projections are not the same size!\n')
        return;
    end
end

% initialize skip counter for unmatched projections
skipcount = 0;

for k = 1:size(proj_h,3)
    if (~mod(k,50))
        fprintf("Decomposing: projection %d/%d\n", k, size(proj_h,3));
    end

    L_array = proj_interp(proj_l, angles_l, angles_h(k), model, tolerance);

    if(isnan(L_array))
        fprintf("No matches for projection %d.  Skipping.\n", k);
        skipcount = skipcount+1;
        continue;
    end

    % map pixels from H and L to equivalent thicknesses
    for j = 1:size(proj_h,2)
        for i = 1:size(proj_h,1)
            H = proj_h(i,j,k);
            L = L_array(i,j);
            material1_proj(i,j,k-skipcount) = c(1)*L + c(2)*H + c(3)*(L.^2) + c(4)*L.*H + c(5)*(H.^2) + c(6)*(L.^3) + c(7)*(H.*(L.^2)) + c(8)*(L.*(H.^2)) + c(9)*(H.^3);
            material2_proj(i,j,k-skipcount) = d(1)*L + d(2)*H + d(3)*(L.^2) + d(4)*L.*H + d(5)*(H.^2) + d(6)*(L.^3) + d(7)*(H.*(L.^2)) + d(8)*(L.*(H.^2)) + d(9)*(H.^3);
            angles_out(k-skipcount) = angles_h(k);
        end
    end
end
end

function [scanner,calibration_pair, model, tolerance, material1, material2] = parse_inputs(argin)
opts = {'scanner','kVl', 'kVh', 'interpolation', 'tolerance','material1','material2'};
defaults=ones(length(opts),1);

% check if option has been passed as input
nVarargs = length(argin);
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option is not default, then extract value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end

        if isempty(ind)
            error('TIGRE:DEDecompose:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
        end
        val=argin{jj};
    end

    switch opt
        case 'scanner'
            if default
                warning("You have not given a scanner for DE calibration.\n Dual Energy Calibration is strictly dependant on scanner calibration, so do not use the results here for quantitative analysis.\n Assuming VarianTrueBeam")
                scanner='VarianTrueBeam';
            else
                if ~ischar( val)
                    error('TIGRE:DEDecompose:InvalidInput','Scanner must be a string')
                end
                scanner=val;
            end
        case 'kVl'
            if default
                kVl=80;
            else
                if ~isscalar(val)
                    error('TIGRE:DEDecompose:InvalidInput','KVl must be a scalar')
                end
                kVl=val;
            end
        case 'kVh'
            if default
                kVh=140;
            else
                if ~isscalar(val)
                    error('TIGRE:DEDecompose:InvalidInput','KVh must be a scalar')
                end
                kVh=val;
            end
        case 'interpolation'
             if default
                model='nearest';
            else
                if ~ischar( val)
                    error('TIGRE:DEDecompose:InvalidInput','Interpolation must be a string')
                end
                model=val;
             end
        case 'tolerance'
            if default
                tolerance=1;
            else
                if ~isnumeric(val)
                    error('TIGRE:DEDecompose:InvalidInput','tolerance must be a number')
                end
                tolerance=val;
            end
                case 'material1'
            if default
                warning('Assuming Al for material1')
                material1='Al';
            else
                if ~ischar( val)
                    error('TIGRE:DEDecompose:InvalidInput','material1 must be a string')
                end
                material1=val;
             end
        case 'material2'
            if default
                warning('Assuming PMMA for material2')
                material2='PMMA';
            else
                if ~ischar( val)
                    error('TIGRE:DEDecompose:InvalidInput','material2 must be a string')
                end
                material2=val;
             end
              
        otherwise
            error('TIGRE:DEDecompose:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in DEDecompose()']);
    end
end
calibration_pair=['V_', int2str(kVl), '_' , int2str(kVh)];
end
