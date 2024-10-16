function [proj,geo]=flatten_detector(proj,geo)
    % This function creates projections on a flat panel detector from a curved
    
    % Check if 'oversample' is provided, otherwise set to default value of 1
    if nargin < 3
        oversample = 1; % Default value for oversample
    end
    
    % Assume we know the arc-length (this has to be input or computed somehow).
    % We assume this is to the end of the pixel, not to the center of the last pixel
    arclength = geo.sDetector(1);
    pixel_arclength = arclength / size(proj, 2);
    
    
    % Calculate the size of the detector projected from the arc onto the plane
    % This is slightly different for odd and even detectors
    if mod(size(proj,2),2)
        % Odd number of detectors - one on the centreline
        detector_size = tan(pixel_arclength) * geo.DSD;
    else
        detector_size = tan(pixel_arclength / 2) * geo.DSD * 2;
    end
    detector_size = detector_size / oversample;
    total_detector_size = tan(arclength / 2) * geo.DSD * 2;

    % Calculate the equal distance detector positions on the plane detector
    if mod(orig_num_detectors,2) 
        % Odd number of detectors - one on the centreline
        detectors = detector_size:detector_size:(total_detector_size / 2);
        detectors = [-fliplr(detectors), 0, detectors];
    else
        detectors = (detector_size / 2):detector_size:(total_detector_size / 2);
        detectors = [-fliplr(detectors), detectors];
    end
    num_detector_columns = length(detectors);

    % Calculate the angle from the source to detector positions on the plane detector
    angles = atan2(detectors, geo.DSD);

    % Normailse to give this angle relative to the projection detector
    normalised_angles = angles / pixel_arclength + orig_num_detectors / 2;

    % Interpolate for projected detectors
    flattened_proj = zeros(size(proj), 'single');
    for i = 1:num_angles
        for r = 1:num_detector_rows
            flattened_proj(i, r, :) = interp1(0:(orig_num_detectors-1), squeeze(proj(i, r, :)), normalised_angles);
        end
    end

    % Update geometry
    geo.nDetector = [num_detector_rows, num_detector_columns];
    geo.dDetector = [geo.dVoxel(1), detector_size];
    geo.sDetector = geo.nDetector .* geo.dDetector;
   
    
end