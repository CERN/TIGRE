function geo=spotDriftGeo(geo,drift,angles)
% SPOTDRIFTGEO geometry for a source focal spot that drifts during the scan
%
%   geo=spotDriftGeo(geo,drift,angles) describes a source that moves away
%   from its nominal position. drift is 2xN or 3xN, one column per
%   projection, or a single column for a constant offset. Rows, in mm, in the
%   frame that rotates with the scanner: drift along u (same direction as
%   geo.offDetector(1,:)), along v (same direction as geo.offDetector(2,:),
%   +z for a standard scan) and, optionally, along the central ray towards
%   the detector. angles are the projection angles, 1xN or 3xN, as given to
%   Ax, Atb and the algorithms.
%
%   The drift is described exactly with fields TIGRE already has: the source
%   and the detector are moved back together, so the source returns to its
%   nominal place, geo.offDetector takes the opposite shift, and
%   geo.offOrigin takes it too, rotated with each angle. A drift along the
%   central ray changes geo.DSO and geo.DSD. All of these become
%   per-projection values.
%
%   See also tiltedAxisGeo
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
assert(~(isfield(geo,'mode') && strcmp(geo.mode,'parallel')),'TIGRE:spotDriftGeo:InvalidInput','spotDriftGeo needs a cone-beam geometry; a parallel beam has no source position.');
if size(angles,1)==1
    angles=[angles;zeros(2,size(angles,2))];
end
assert(size(angles,1)==3,'TIGRE:spotDriftGeo:InvalidInput','angles should be 1xN or 3xN.');
n=size(angles,2);
assert(any(size(drift,1)==[2 3]) && any(size(drift,2)==[1 n]),'TIGRE:spotDriftGeo:InvalidInput','drift should be 2xN or 3xN, or a single column, with N the number of angles.');
drift=repmat(drift,1,n/size(drift,2));
dU=drift(1,:); dV=drift(2,:);
if size(drift,1)==3
    dW=drift(3,:);
else
    dW=zeros(1,n);
end

% Moving source and detector by -drift brings the source back; the image moves with them.
% The rotating frame's u and v axes in world coordinates (columns 2 and 3 of the Euler ZYZ
% rotation the projectors use).
a=angles(1,:); t=angles(2,:); p=angles(3,:);
axU=[-cos(a).*cos(t).*sin(p)-sin(a).*cos(p); -sin(a).*cos(t).*sin(p)+cos(a).*cos(p); sin(t).*sin(p)];
axV=[cos(a).*sin(t); sin(a).*sin(t); cos(t)];
if ~isfield(geo,'offOrigin') || isempty(geo.offOrigin)
    geo.offOrigin=[0;0;0];
end
if ~isfield(geo,'offDetector') || isempty(geo.offDetector)
    geo.offDetector=[0;0];
end
geo.offOrigin=perAngle(geo.offOrigin,n)-(dU.*axU+dV.*axV);
geo.offDetector=perAngle(geo.offDetector,n)-[dU;dV];
geo.DSO=perAngle(geo.DSO,n)-dW;
geo.DSD=perAngle(geo.DSD,n)-dW;
end

function v=perAngle(v,n)
if size(v,2)==1
    v=repmat(v,1,n);
end
assert(size(v,2)==n,'TIGRE:spotDriftGeo:InvalidInput','geometry fields should have one column, or one per angle.');
end
