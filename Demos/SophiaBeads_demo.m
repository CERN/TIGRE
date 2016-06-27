% This demo Will show how to reconstruct real data. 
%
% We will reconstruct the so-called SophiaBeads dataset, where a bunch of
% "balls" are put in a Nikon ?-CT machine.
%
% The data can be downloaded at:
% https://zenodo.org/record/16474
%
% And the codes can be downloaded at:
% http://zenodo.cern.ch/record/16539
%
% The SophiaBeads codes and dataset ARE NOT part of TIGRE toolbox, and are
% not licensed under the BSD license, be aware.
%
% Additionally, if you use the dataset or the codes, they have each of them
% their own bilbiography, that you need to reference, as per their own
% requirements. 

% cite as:

% citeme('SophiaBeads data')
% citeme('SophiaBeads code')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------

%% How to run this code:
%
% 1.- Download both the code and the data (any size you prefer). 
%
% 2.- Open SophiaBeads.m in the codes and replace pathname and filename (in
%     lines 13-14 by the data path and name you chosen
%
% 3.- In line 27 add the followinh code:     break; 
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
% very big size and most ocmputers will not be able to run it. We redefine
% it in the next line, but feel free to use the original
geo.nVoxel=[200 200 200].';
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


 data=permute(data,[2 1 3]);

%% Now use any algorithms you want

% example:
beadsSART=SART(data,geo,angles,15);
beadsCGLS=CGLS(data,geo,angles,15);

plotImg([beadsSART beadsCGLS],'Dim','Z','Step',1)
