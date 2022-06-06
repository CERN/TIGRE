function nmask = SC_GroupMask(thickness, ngroup, nbounds)
%% Generate n-group masks for estimated thickness (2D)
%
% SYNOPSIS: nmask = GroupMask(thickness, nbounds)
%
% INPUT nbounds: thickness based n-group bounds
%		thickness: estimated thickness 2D matrix
%
% OUTPUT nmask(u, v, ngroup)
%
% REMARKS
%
% created with MATLAB ver.: 9.8.0.1538580 (R2020a) Update 6 on Microsoft Windows 10 Pro Version 10.0 (Build 18363)
%
% created by: Yi Du
% DATE: 18-May-2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% thickness based n-group
nmask = zeros(size(thickness, 1), size(thickness, 2), ngroup);

%% n-group mask
for ii = 1:ngroup-1
    nmask(:, :, ii) = (thickness>nbounds(ii)).*(thickness<nbounds(ii+1));
end

nmask(:,:,ngroup) = thickness>nbounds(ngroup);

end