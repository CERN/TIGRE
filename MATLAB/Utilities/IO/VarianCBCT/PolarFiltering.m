function img3D_f = PolarFiltering(img3D)
% Polar Filtering in image domain for ring correction
% Method: dependent on third-party tool, polar2cart
%       ring artifacts is first to be identified, and only one ring can be
%       corrected in each run.
% Input:
%       img3D: image matrix
%       geo: geometry
% Output:
%       img: filtered projection matrix
% Date: 2022-02-04
% Author: Yi Du, yi.du@hotmail.com


imgSize = size(img3D, 1);
img3D_f = img3D;

%% slice by slice
for ii = 1:size(img3D,3)    
    % Cartisan to Polar
    imP = ImToPolar(img3D(:,:,ii), 0, 1, imgSize, imgSize);
    sum_over_col = sum(imP,2);
    
    % find the indices of max and min values
    [~, max_idx] = max(diff(sum_over_col), [], 'linear');    
    [~, min_idx] = min(diff(sum_over_col), [], 'linear');
    
    % median filtering
    if(abs(max_idx - min_idx)>10)
        continue;
    else
        if(max_idx > min_idx)
            imP(min_idx:max_idx,:) = ordfilt2(imP(min_idx:max_idx,:), 5, ones(5,1));
        else
            imP(max_idx:min_idx,:) = ordfilt2(imP(max_idx:min_idx,:), 5, ones(5,1));
        end
    end
    
    img3D_f(:,:,ii) = PolarToIm(imP, 0, 1, imgSize, imgSize);
    
end

end

