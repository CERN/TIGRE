function [ img ] = Atb( projections,geo,angles,varargin )
%ATB CUDA backprojection operator

% Lets make 100% sure that data is correct
%% OPtionals

ptype='FDK';
if nargin > 3
   assert(any(strcmpi(varargin{1},{'FDK','matched'})),'TIGRE:Atb:InvalidInput','Projection type not understood (4th input).');
   ptype=varargin{1};
end


%% image
assert(isa(projections,'single'),'TIGRE:Atb:InvalidInput','Image shoudl be single type');
assert(isreal(projections),'TIGRE:Atb:InvalidInput','Image shoudl be real (non-complex)');
assert(size(projections,2)>1,'TIGRE:Atb:InvalidInput', 'Projections shoudl be 2D'); %TODO: needed? 
assert(size(projections,3)==size(angles,2),'TIGRE:Atb:InvalidInput', 'Number of projections shoudl match number of angles.'); 
%% Angles
assert(isreal(angles),'TIGRE:Atb:InvalidInput','Angles shoudl be real (non-complex)');
assert(size(angles,1)==1 | size(angles,1)==3 ,'TIGRE:Atb:InvalidInput','Angles shoudl be of size 1xN or 3xN');
angles=double(angles); %in case they were single.

%% geometry
geo=checkGeo(geo,angles);
assert(isequal([size(projections,2) size(projections,1)],geo.nDetector.'),'TIGRE:checkGeo:BadGeometry','nVoxel does not match with provided image size');

%% Thats it, lets call the mex fucntion

img=Atb_mex(projections,geo,angles,ptype);



end

