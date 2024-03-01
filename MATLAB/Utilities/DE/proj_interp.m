function proj_out = proj_interp(proj_in, angles_in, target_angle, model, tolerance)
%PROJ_INTERP Performs sinogram interpolation of the input projections
% according to the specified interpolation model
%   proj_out is the interpolated projection at target_angle
%
%   proj_in is the list of projections to be interpolated, given
%   at angles listed in angles_in 
%
%   model specifies the interpolation model to be used.
%       'linear' uses a linear interpolation model
%       'nearest' uses a nearest-neighbor interpolation model
%       'skip' is the same as 'nearest' but will skip projections without
%       sufficiently close matches. Requires setting 'tolerance.'
%   
%   Coded by: Andrew Keeler
%   Contact: akeeler@luc.edu
%   Date: February 28, 2024

[delta, idx, idx_next] = angle_align(target_angle, angles_in);

switch model
    case 'linear'
        dy = proj_in(:,:,idx) - proj_in(:,:,idx_next);
        proj_out = proj_in(:,:,idx) + dy*delta;
    case 'nearest'
        proj_out = squeeze(proj_in(:,:,idx));
    case 'skip'
        if delta > tolerance
            proj_out = NaN;
            return;
        end
        proj_out = squeeze(proj_in(:,:,idx));
end

end

function [weight, idx, idx_next] = angle_align(target, angle_array)

angle_diff = angle_array - target;
last = length(angle_diff);
[~, idx] = min(abs(angle_diff));

if (angle_diff(idx) == 0)
    idx_next = idx;
    weight = 1;
    return
end

if (idx == 1 && sign(angle_diff(1)) == sign(angle_diff(2)))
    idx_next = idx;
    weight = 1;
    return
elseif (idx == 1 && sign(angle_diff(1)) ~= sign(angle_diff(2)))
    idx_next = 2;
    weight = angle_diff(idx)/(abs(angle_array(idx) - angle_array(idx_next)));
    return
end



if (idx == last && sign(angle_diff(last)) == sign(angle_diff(last-1)))
    idx_next = idx;
    weight = 1;
    return
end

if (sign(angle_diff(idx)) ~= sign(angle_diff(idx-1)))
    idx_next = idx - 1;
    weight = angle_diff(idx)/(abs(angle_array(idx) - angle_array(idx_next)));
else
    idx_next = idx + 1;
    weight = angle_diff(idx)/(abs(angle_array(idx) - angle_array(idx_next)));
end
end