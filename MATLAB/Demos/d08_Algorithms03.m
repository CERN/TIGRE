%% DEMO 8:  Algorithms 03. Krylov subspace
%
%
% In this demo the usage of the the Krylov subspace family is explained.
% This family of algorithms iterates trhough the eigenvectors of the
% residual (Ax-b) of the problem in descending order, achieving increased
% convergence rates comparing to SART family. 
% 
% In cases where the data is good quality, SART type families tend to reach
% to a better image, but when the data gets very big, or has bad quality,
% CGLS is a good and fast algorithm
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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------
%% Initialize
clear;
close all;

%% Define Geometry
geo=defaultGeometry('nVoxel',[512,512,512]','nDetector',[512,512]);                     

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi,100);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);

%% Usage CGLS
%
%
%  CGLS has the common 4 inputs for iterative algorithms in TIGRE:
%
%  Projections, geometry, angles, and number of iterations 
%
% Additionally it contains optional initialization tehcniques, but we
% reccomend not using them. CGLS is already quite fast and using them may
% lead to divergence.
% The options are:
%  'Init'    Describes diferent initialization techniques.
%             �  'none'     : Initializes the image to zeros (default)
%             �  'FDK'      : intializes image to FDK reconstrucition
%             �  'multigrid': Initializes image by solving the problem in
%                            small scale and increasing it when relative
%                            convergence is reached.
%             �  'image'    : Initialization using a user specified
%                            image. Not recomended unless you really
%                            know what you are doing.
%  'InitImg'    an image for the 'image' initialization. Avoid.
 
% use CGLS
[imgCGLS, errL2CGLS]=CGLS(noise_projections,geo,angles,60);
% SIRT for comparison.
[imgSIRT,errL2SIRT]=SIRT(noise_projections,geo,angles,60);

%% plot results
%
% We can see that CGLS gets to the same L2 error in less amount of
% iterations.

plot([errL2SIRT;[errL2CGLS nan(1,length(errL2SIRT)-length(errL2CGLS))]]');
title('L2 error')
legend('SIRT','CGLS')

% plot images
plotImg([imgCGLS imgSIRT],'Dim','Z','Step',2)
%plot errors
plotImg(abs([head-imgCGLS head-imgSIRT]),'Dim','Z','Slice',64)
