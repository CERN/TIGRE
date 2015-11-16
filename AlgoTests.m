%% ALGOTESTS
% This function tests different algorithms of image recosntruction
% using the CUDA CBCT code.
% 
% Add new test of functionallity here, for future reference.


%% Set up enviroment.

clear;
clc;
close all;

% Perceptually uniform colormaps
addpath('.\Colormaps');
% Add third party tools from FEX
addpath('.\Third_party_tools\3DSheppLogan'); % 3D shepp-Logan
addpath('.\Third_party_tools\arro3d'); % 3D shepp-Logan
% Add tolbox folders
addpath('.\Algorithms');
addpath('.\Utilities');
%% Create geometry

% Load a default geometry (look inside if curious)
GeometrySetting;

% Set accuracy (1/amount of samples per voxel)
Geometry.accuracy=0.5;

% Set projection angles
init=0;          % start angle
step=0.5;        % step
finish=360-step; % end angle

alpha=[init:step:finish]*pi/180;

%% Plot the geometry

% plotgeometry(Geometry,30)

%% Define the test image.




%% Generate some data


