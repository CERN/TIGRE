 function centre=computeCOR(data,geo,angles,slice)
% Code modified from SophiaBeads
% (https://github.com/Sophilyplum/sophiabeads-datasets/blob/master/tools/centre_geom.m)
% Reference:
% T. Liu - "Direct central ray determination in computed microtomography",
% Optical Engineering, April 2009.

% Get central slice
if nargin<4
    slice=floor(size(data,1)/2)+1;
end
data=squeeze(data(slice,:,:));

if size(angles,1)==1
   angles=angles'; 
end

% Set up coordinate grids for testing the fit to data
[angle_grid, det_grid] = meshgrid(angles',linspace(-geo.sDetector(1)/2+geo.dDetector(1)/2,+geo.sDetector(1)/2-geo.dDetector(1)/2,geo.nDetector(1))');
% Wrap grids in the angular direction to avoid value out of range errors
angle_grid = [(angle_grid - 2*pi) angle_grid (angle_grid + 2*pi)];
det_grid = repmat(det_grid,1,3);
test_data = double(repmat(data,1,3));

% Start search using midpoint at zero
midpoint = 0;
% Vector of precision values to search at
precision = [1 0.1 0.01 0.001 0.0001 0.00001 0.000001];

for i = 1:length(precision)
    
    COR = (midpoint - 10*precision(i)):precision(i):(midpoint + 10*precision(i)); % values for centre of rotation
    M = zeros(length(COR),1);
    
    for j = 1:length(COR)
        
        gamma = atan(linspace(-geo.sDetector(1)/2+geo.dDetector(1)/2,+geo.sDetector(1)/2-geo.dDetector(1)/2,geo.nDetector(1))' / geo.DSD);   % angle of each ray relative to theoretical central ray    
        gamma_c = atan(COR(j) / geo.DSD); % angle of assumed centre of rotation to central ray
        gamma_i = gamma - gamma_c;
        beta = 2 * gamma_i + pi;
        
        s2 = geo.DSD * tan(2 * gamma_c - gamma);
        s2 = repmat(s2, 1, size(angles,2));
        
        angles_aux = repmat(angles', geo.nDetector(1), 1) + repmat(beta, 1, size(angles,2));
        test = interp2(angle_grid, det_grid, test_data, angles_aux, s2, 'linear', 0);
        
        nonzero = find(test > 0);
        % We want the number of non-zero values for the average, not the sum of their positions
        M(j) = sum((test(nonzero) - data(nonzero)).^2)*(1/length(nonzero));
    end
    
    [~, indM] = min(M);   % minimum value and index
    midpoint = COR(indM);   % set midpoint for next search
end

% transform centre to required value
centre =- midpoint * geo.DSO / geo.DSD;

