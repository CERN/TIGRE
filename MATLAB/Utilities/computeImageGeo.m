function geo=computeImageGeo(geo)
% COMPUTEIMAGEGEO image size and voxel size that suit the detector and the scan
%
%   geo=computeImageGeo(geo) sets geo.nVoxel, geo.dVoxel and geo.sVoxel. The
%   voxels get the size of a detector pixel seen at the rotation axis,
%   dDetector*DSO/DSD (dDetector for parallel beam), and there are enough of
%   them to cover the field of view: the circle that every projection sees in
%   the axial plane (larger with an offset, half-fan, detector), and the
%   detector's height at the axis. Useful for real data, where the image size
%   is not given. Per-projection DSO and DSD are averaged, and the largest
%   offDetector is used.
%
%   See also changeResolution
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
dU=geo.dDetector(1); dV=geo.dDetector(2);
offU=0; offV=0;
if isfield(geo,'offDetector') && ~isempty(geo.offDetector)
    offU=max(abs(geo.offDetector(1,:)));
    offV=max(abs(geo.offDetector(2,:)));
end
uMax=geo.nDetector(1)*dU/2+offU; % furthest detector column from the central ray, mm
vMax=geo.nDetector(2)*dV/2+offV;
if isfield(geo,'mode') && strcmp(geo.mode,'parallel')
    mag=1;
    radius=uMax;
else
    DSO=mean(geo.DSO); DSD=mean(geo.DSD);
    mag=DSD/DSO;
    radius=DSO*uMax/hypot(DSD,uMax); % distance of the outermost ray from the axis
end
dXY=dU/mag;
dZ=dV/mag;
nXY=ceil(2*radius/dXY);
nZ=ceil(2*vMax/mag/dZ);
geo.nVoxel=[nXY;nXY;nZ];
geo.dVoxel=[dXY;dXY;dZ];
geo.sVoxel=geo.nVoxel.*geo.dVoxel;
end
