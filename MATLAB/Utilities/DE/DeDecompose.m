function [Al_thickness,PMMA_thickness, angles_out] = DeDecompose(proj_l, angles_l, proj_h, angles_h, varargin)
%DEDECOMPSE Maps CBCT projections taken at two kVp settings to equivalent
%thicknesses of Al and PMMA using precalibrated transfer functions.
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

if ~ isfile('./de_cal.mat')
    warning('Using DEMO Dual-energy transfer calibration. We recommend creating your own and storing it in TIGRE/Common/data/de_cal.mat')
    load('./demo_de_cal.mat', 'de_cal');
else
    load('./de_cal.mat', 'de_cal');
end

[calibration, model, tolerance] = parse_inputs(varargin{:});

tmp = getfield(de_cal, calibration);

c = tmp.Al;
d = tmp.PMMA;

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
    
    if(L_array == NaN)
        fprintf("No matches for projection %d.  Skipping.\n", k);
        skipcount = skipcount+1;
        continue;
    end

    % map pixels from H and L to equivalent thicknesses
    for j = 1:size(proj_h,2)
        for i = 1:size(proj_h,1)
            H = proj_h(i,j,k);  
            L = L_array(i,j);
            Al_thickness(i,j,k-skipcount) = c(1)*L + c(2)*H + c(3)*(L.^2) + c(4)*L.*H + c(5)*(H.^2) + c(6)*(L.^3) + c(7)*(H.*(L.^2)) + c(8)*(L.*(H.^2)) + c(9)*(H.^3);
            PMMA_thickness(i,j,k-skipcount) = d(1)*L + d(2)*H + d(3)*(L.^2) + d(4)*L.*H + d(5)*(H.^2) + d(6)*(L.^3) + d(7)*(H.*(L.^2)) + d(8)*(L.*(H.^2)) + d(9)*(H.^3);
            angles_out(k-skipcount) = angles_h(k);
        end
    end
end
end

function [calibration, model, tolerance] = parse_inputs(varargin)
    opts = {'kVl', 'kVh', 'model', 'tolerance'};
    args = {80, 140, 'nearest', 1};

    % Check inputs
    nVarargs = length(varargin);
    if mod(nVarargs,2)
        error('TIGRE:FDK:InvalidInput','Invalid number of inputs')
    end

    % check if option has been passed as input, assign new value if it has
    for ii=1:2:nVarargs
        ind=find(ismember(opts,varargin{ii}));
        if ~isempty(ind)
            args{ind}=varargin{ii+1};
        end
    end

    for jj = 1:length(opts)
        opt = opts{jj};
        switch opt
            case 'kVl'
                kVl = args{jj};
            case 'kVh'
                    kVh = args{jj};
            case 'model'
                    model = args{jj};
            case 'tolerance'
                    tolerance = args{jj};
        end
    end
    calibration = "V_"+num2str(kVl)+"_"+num2str(kVh);
end