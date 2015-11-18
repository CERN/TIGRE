%% ALGOTESTS
% This function tests different algorithms of image recosntruction
% using the CUDA CBCT code.
% 
% Add new test of functionallity here, for future reference.


%% Set up enviroment.

clear;
clc;
close all;


% Add tolbox folders
addpath('.\Algorithms');
addpath('.\Utilities');
addpath('.\Test_data');
addpath('.\Mex_files');

% Perceptually uniform colormaps
addpath('.\Colormaps');
% Add third party tools from FEX
addpath('.\Third_party_tools\arrow3d'); % 3D shepp-Logan
%% Create geometry

% Load a default geometry (look inside if curious)
GeometrySetting;

% Set accuracy (1/amount of samples per voxel)
Geometry.accuracy=0.5;

% Set projection angles
init=0;          % start angle
step=1;        % step
finish=360-step; % end angle

alpha=[init:step:finish]*pi/180;

%% Plot the geometry

plotgeometry(Geometry,30)

%% Define the test image.
 
% 3D Shepp-Logan
% img=phantom3dAniso('Modified Shepp-Logan',Geometry.nVoxel);

% thorax digital phantom
img=thoraxPhantom(Geometry.nVoxel);

%% Plot the image

plotImg(img,'Step',1,'Dim',3);

%% Generate some data

data=Ax(img,Geometry,alpha);
% Add noise to data;

maxval=max(data(:));
data_noise=imnoise(data/maxval,'poisson');
data_noise=data_noise.*maxval;
data=data_noise;
clear data_noise
%% Plot the data

plotProj(data,alpha);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIRT 
[resSIRT, errSIRT]=SIRT_CBCT(data,Geometry,alpha,15);
%% SART

% Due to memory management stuff, SART is about 5 times slower
[resSART, errSART]=SART_CBCT(data,Geometry,alpha,15);

%% OS-SART
block_size=10;
[resOSSART, errOSSART]=SART_CBCT(data,Geometry,alpha,10,15);


