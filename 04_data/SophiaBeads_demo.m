% Real micro-CT data: SophiaBeads
% ===============================
%
% Reconstruct the SophiaBeads dataset, a tube of glass beads.
%
% Source: `MATLAB/Demos/SophiaBeads_demo.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/SophiaBeads_demo.m>`__
%
% This demo is only available in MATLAB, so it is shown without output.

%% How to run this code

%
% 1.- Download both the code and the data (any size you prefer).
%
% 2.- Open SophiaBeads.m in the codes and replace pathname and filename (in
%     lines 13-14 by the data path and name you chosen
%
% 3.- In line 27 add the following code:     break;
%
% 4.- Run SophiaBeads.m
%
% 5.- Without clearing the variables, come to this code and run it.
%
%

% from sophiabeads to TIGRE
geo.DSD = geom.d_sd;
geo.DSO = abs(geom.source.x);

% Detector parameters
geo.nDetector=[geom.dets.ny;geom.dets.nz];
geo.dDetector=[mean(diff(geom.dets.y)); mean(diff(geom.dets.z))];
geo.sDetector=geo.nDetector.*geo.dDetector;

% Image parameters
geo.nVoxel=geom.voxels.';  % this is the datasheets standard size, but its
% very big size and most computers will not be able to run it. We redefine
% it in the next line, but feel free to use the original
geo.nVoxel=[500 500 200].';
geo.sVoxel=[geom.voxel_size.*geom.voxels].';
geo.dVoxel=geo.sVoxel./geo.nVoxel;

% Offsets
geo.offOrigin=[0,0,0].';
geo.offDetector=[0; 0];

% Auxiliary
geo.accuracy=0.5;
geo.mode='cone';
geo.COR=geom.source.y;

% Angles
angles=-geom.angles.';

data=reshape(data,1564,1564,size(angles,2));
data=permute(data,[2 1 3]);

%% Now use any algorithms you want

% example:
% beadsSART=SART(data,geo,angles,15);
beadsCGLS=CGLS(data,geo,angles,15);
ImgASD_POCS_20=OS_ASD_POCS(data,geo,angles,50,...
'TViter',2,'maxL2err',0.5,'alpha',0.002,'alpha_red',0.2,... % these are very important
'lambda',1,'lambda_red',0.99,'Ratio',1,'Verbose',true);
% plotImg([beadsSART beadsCGLS],'Dim','Z','Step',1)
