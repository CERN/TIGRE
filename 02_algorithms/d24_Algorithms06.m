% Prior images with PICCS
% =======================
%
% Use a prior image of the same object, such as a planning CT.
%
% Source: `MATLAB/Demos/d24_Algorithms06.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d24_Algorithms06.m>`__
%
% This demo is only available in MATLAB, so it is shown without output.

%% Initialize

clear;
close all;

%% Define Geometry

geo=defaultGeometry('nVoxel',[128;128;128]);

%% Load data and generate projections

% see previous demo for explanation
angles=linspace(0,2*pi-2*pi/30,30);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);

%% Lets create a OS-SART test for comparison

[imgOSSART,errL2OSSART]=OS_SART(noise_projections,geo,angles,60);

%% ASD\_POCS (TV without prior)

epsilon=errL2OSSART(end);
alpha=0.002;
ng=25;

imgASDPOCS=ASD_POCS(noise_projections,geo,angles,50,'TViter',ng,'maxL2err',epsilon,'alpha',alpha);

%% PICCS (ASD\_POCS, but with an image prior)

imgPICCS=PICCS(noise_projections,geo,angles,50,head,'TViter',ng,'maxL2err',epsilon,'alpha',alpha);

%% Plot results

errOSSART=abs(imgOSSART-head);
errASDPOCS=abs(imgASDPOCS-head);
errPICCS=abs(imgPICCS-head);

plotImg( [imgOSSART imgASDPOCS imgPICCS],'Dim','z','slice',34);
plotImg( [errOSSART errASDPOCS errPICCS],'Dim','z','slice',34,'clims',[0,0.1]);
