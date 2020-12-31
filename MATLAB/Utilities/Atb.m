function [ img ] = Atb( projections,geo,angles,varargin )
%ATB CUDA backprojection operator
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
% Codes:              https://github.com/CERN/TIGRE/blob/master/LICENSE
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
% Lets make 100% sure that data is correct
%% OPtionals

ptype='FDK';
if nargin > 3
   assert(any(strcmpi(varargin{1},{'FDK','matched'})),'TIGRE:Atb:InvalidInput','Projection type not understood (4th input).');
   ptype=varargin{1};
end


%% image
assert(isa(projections,'single'),'TIGRE:Atb:InvalidInput','Image should be single type');
assert(isreal(projections),'TIGRE:Atb:InvalidInput','Image should be real (non-complex)');
assert(size(projections,2)>1,'TIGRE:Atb:InvalidInput', 'Projections should be 2D'); %TODO: needed? 
assert(size(projections,3)==size(angles,2),'TIGRE:Atb:InvalidInput', 'Number of projections should match number of angles.'); 
%% Angles
assert(isreal(angles),'TIGRE:Atb:InvalidInput','Angles should be real (non-complex)');
assert(size(angles,1)==1 | size(angles,1)==3 ,'TIGRE:Atb:InvalidInput','Angles should be of size 1xN or 3xN');
angles=double(angles); %in case they were single.
if size(angles,1)==1
   angles=repmat(angles,[3 1]);
   angles(2,:)=0;
   angles(3,:)=0;
end
%% geometry
geo=checkGeo(geo,angles);
assert(isequal([size(projections,2) size(projections,1)],geo.nDetector.'),'TIGRE:checkGeo:BadGeometry','nVoxel does not match with provided image size');

%% Thats it, lets call the mex fucntion

img=Atb_mex(projections,geo,angles,ptype);



end

