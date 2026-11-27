function geo=tiltedAxisGeo(geo,corTop,corBottom,height)
% TILTEDAXISGEO geometry for a rotation axis tilted within the detector plane
%
%   geo=tiltedAxisGeo(geo,corTop,corBottom) describes a rotation axis that is
%   not parallel to the detector columns: its centre of rotation is corTop at
%   the top of the image volume (+z) and corBottom at the bottom (-z), for
%   example two values of computeCOR at two detector rows. The top is the end
%   seen by the last detector rows, proj(end,:,:); row 1 is the bottom. Both
%   values are in mm at the rotation axis, with the sign that geo.COR uses.
%
%   geo=tiltedAxisGeo(geo,corTop,corBottom,height) gives the distance between
%   the two measurements along the rotation axis, in mm. Default: geo.sVoxel(3).
%
%   The tilt is described exactly with fields TIGRE already has: geo.COR is
%   replaced, the tilt angle is added to the in-plane detector rotation
%   geo.rotDetector(1,:), and a small shift is added to geo.offOrigin(3,:).
%   The reconstructed image is aligned with the rotation axis.
%
%   Iterative algorithms are exact. FDK filters along the detector rows, so
%   it is approximate for a tilted axis; for large tilts, rotate the
%   projections by the tilt angle instead and set only geo.COR.
%
%   See also computeCOR
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
if nargin<4 || isempty(height)
    height=geo.sVoxel(3);
end
assert(isscalar(corTop) && isscalar(corBottom) && isscalar(height),'TIGRE:tiltedAxisGeo:InvalidInput','corTop, corBottom and height should be scalars.')
assert(height>0,'TIGRE:tiltedAxisGeo:InvalidInput','height should be positive.')
s=(corBottom-corTop)/height;
assert(abs(s)<1,'TIGRE:tiltedAxisGeo:InvalidInput','corTop and corBottom differ by more than height.')

phi=asin(s);                 % tilt of the axis within the detector plane
corMid=(corTop+corBottom)/2; % centre of rotation half way between the two

geo.COR=corMid*cos(phi);
if ~isfield(geo,'rotDetector')
    geo.rotDetector=[0;0;0];
end
geo.rotDetector(1,:)=geo.rotDetector(1,:)+phi;
if ~isfield(geo,'offOrigin')
    geo.offOrigin=[0;0;0];
end
geo.offOrigin(3,:)=geo.offOrigin(3,:)-corMid*sin(phi);
end
