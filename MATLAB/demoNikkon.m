
clear;clc;


filepath='put the path to the data!'; % e.g. C:\Mydata\phantom1     


% Just add the path. This assumes that a .xtekct adn a .ang file is inside
% that folder
[geo,angles,whiteLevel]=readXtekctGeometry(filepath);

% This function will load and preprocess the data. 
% call it as 
% [loadeddata,resampled_angles]=loadNikkonData(path,all_the_angles,number_of_desired_samples,whiteLevel)
[proj,angles]=loadNikkonData(filepath,angles,20,whiteLevel);

%% plot the loaded data
% plotProj(proj,angles);

%% make the image it a bit smaller (2000x2000x2000 wont fit in memory)
geo.nVoxel=geo.nVoxel/10;
geo.dVoxel=geo.sVoxel./geo.nVoxel;

%% test a random algorithm

img=CGLS(proj,geo,angles,10);

% plot a slice in Z dimension, the middle one (number 100).
plotImg(img,'Dim','z','slice',100)
