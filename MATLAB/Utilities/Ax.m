function [ projections ] = Ax(img, geo, angles, varargin )
%AX computes projections for images and geometry information

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

ptype='ray-voxel';
if nargin > 3
   assert(any(strcmpi(varargin{1},{'ray-voxel','interpolated'})),'TIGRE:Ax:InvalidInput','Projection type not understood (4th input).');
   ptype=varargin{1};
end


%% image
assert(isa(img,'single'),'TIGRE:Ax:InvalidInput','Image shoudl be single type');
assert(isreal(img),'TIGRE:Ax:InvalidInput','Image shoudl be real (non-complex)');
% assert(size(img,3)>1,'TIGRE:Ax:InvalidInput', 'image shoudl be 3D'); %TODO: needed? 
%% Angles
assert(isreal(angles),'TIGRE:Ax:InvalidInput','Angles shoudl be real (non-complex)');
assert(size(angles,1)==1 | size(angles,1)==3 ,'TIGRE:Ax:InvalidInput','Angles shoudl be of size 1xN or 3xN');
angles=double(angles); %in case they were single.
if size(angles,1)==1
   angles=repmat(angles,[3 1]);
   angles(2,:)=0;
   angles(3,:)=0;
end
%% geometry
geo=checkGeo(geo,angles);
assert(isequal([size(img,1) size(img,2) size(img,3)],squeeze(geo.nVoxel.')),'TIGRE:Ax:BadGeometry','nVoxel does not match with provided image size');

%% Temporary (?)

s= sum(abs(angles),2);
if (s(2)~=0 || s(3)~=0) && strcmp(ptype,'ray-voxel') && strcmp(geo.mode,'parallel')
    warning(['''ray-voxel'' Not supported for parallel beam arbitrary axis rotation, changed to ''interpolated''.',char(10),'Ignore this message if you are not purposedly triying enforce its use.']);
    ptype='interpolated';
end


%% Thats it, lets call the mex fucntion



projections=Ax_mex(img,geo,angles,ptype);

