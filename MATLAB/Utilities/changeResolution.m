function [geo,proj]=changeResolution(geo,factor,proj)
% CHANGERESOLUTION coarser version of a geometry, and optionally of its projections
%
%   geo=changeResolution(geo,factor) keeps the image's physical size,
%   geo.sVoxel, with factor times fewer voxels along each axis.
%
%   [geo,proj]=changeResolution(geo,factor,proj) also bins the projections by
%   the same factor (the mean of each factor x factor block of pixels) and
%   changes the detector to match. Pixels that do not fill a whole block, at
%   the end of each detector row or column, are dropped, and offDetector moves
%   so that the binned pixels stay where they were. Useful to try parameters
%   on a small problem first.
%
%   See also computeImageGeo
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
assert(isscalar(factor) && factor>=1 && factor==round(factor),'TIGRE:changeResolution:InvalidInput','factor should be an integer, 1 or more.');
geo.nVoxel=max(round(geo.nVoxel/factor),1);
geo.dVoxel=geo.sVoxel./geo.nVoxel;
if nargin<3
    return
end
nU=geo.nDetector(1); nV=geo.nDetector(2);
assert(size(proj,1)==nV && size(proj,2)==nU,'TIGRE:changeResolution:InvalidInput','proj should be of size [nDetector(2) nDetector(1) nangles].');
nUnew=floor(nU/factor); nVnew=floor(nV/factor);
assert(nUnew>=1 && nVnew>=1,'TIGRE:changeResolution:InvalidInput','factor is larger than the detector.');
nAngles=size(proj,3);
kept=proj(1:nVnew*factor,1:nUnew*factor,:);
kept=reshape(kept,factor,nVnew,factor,nUnew,nAngles);
proj=single(reshape(mean(mean(double(kept),1),3),nVnew,nUnew,nAngles));
% The dropped pixels are at the high-index end; the kept ones are centred lower by half of them.
% u grows with the column and v with the row.
shift=[(nU-nUnew*factor)*geo.dDetector(1); (nV-nVnew*factor)*geo.dDetector(2)]/2;
if ~isfield(geo,'offDetector') || isempty(geo.offDetector)
    geo.offDetector=[0;0];
end
geo.offDetector=geo.offDetector-shift;
geo.nDetector=[nUnew;nVnew];
geo.dDetector=geo.dDetector*factor;
geo.sDetector=geo.nDetector.*geo.dDetector;
end
