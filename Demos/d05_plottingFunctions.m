%% Demo 5: How to use plotting functions
%
% This demo will demostrate the options for plotting projection and images
% on TIGRE. The functions have been in previous demos, but in here an
% exaustive explanation and usage of them is given.
%
%
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
%% Initialize

clear;
close all;
%% Define Geometry
% 
% VARIABLE                                   DESCRIPTION                    UNITS
%-------------------------------------------------------------------------------------
geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters
geo.nDetector=[512; 512];					% number of pixels              (px)
geo.dDetector=[0.8; 0.8]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[128;128;128];                   % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi,100);
thorax=thoraxPhantom(geo.nVoxel);
projections=Ax(thorax,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);
%% Reconstruct image using OS-SART and FDK
imgFDK=FDK(noise_projections,geo,angles);
imgOSSART=OS_SART(noise_projections,geo,angles,50);


%% Lets use PlotProj
% 
% plotProj plots the projection data measure on the detector on each angle.
%
% exhaustive list of possible parameters:

% 'Step' : Defines the step size for skippin projections when plotting,
% usefull when there are a big amount of projections. Default is 1
step=2;

% 'Colormap': Defines the colormap used to plot. Default is 'gray'. any
% MATLAB colromap can be used, also any user created Nx3 colorpam can be
% used. Additionally, TIGRE includes perceptually uniform colormaps, where
% ?color == ?data. The perceptuall yuniform options are the default 'gray',
% 'magma', 'viridis', 'inferno' and 'plasma'

colormap=hsv(10);
colormap='gray';
colormap='viridis';
colormap='plasma';
colormap='gray';

% 'Clims': Defines the data limits for the color, usefull when one wants to
% see some specific range. The default uses the minimum and maximum of the
% data.

clims=[0 4];

% 'Savegif': allows to save the plotted figure as an animated gif,
% specified by the given filename.

giffilename='demo5projections.gif';

% 'Slice': allows to plot a single projection .Will overwrite the behaviour
% of 'Step'
slice=5;

% Lets try out. 
plotProj(noise_projections,angles,'Step',step,'Colormap',colormap,'CLims',clims,'Savegif',giffilename); % not using 'Step'

% Remember you can also plot errors, for example the added noise by:

noise=abs(noise_projections-projections); %abs is what we are interested in plotting
plotProj(noise,angles,'CLims',[0,0.05]); 


%% What about plotting the Image? plotImg()
%
% plotImg plots the image slice by slice. 
%
% List of optional parameters: 


% 'Dim': specifies the dimension for plotting. 
%    Dim can be 'X','Y','Z' or 1, 2 ,3 (being the same)
%
dimension='Z'; % or 'z' or 3

% 'Step': step size of the plotting. Useful when images are big or one just
% wants an overview of the result

step=2;
% 'Colormap': Defines the colormap used to plot. Default is 'gray'. any
% MATLAB colromap can be used, also any user created Nx3 colorpam can be
% used. Additionally, TIGRE includes perceptually uniform colormaps, where
% ?color == ?data. The perceptuall yuniform options are the default 'gray',
% 'magma', 'viridis', 'inferno' and 'plasma'

colormap=hsv(20);
colormap='magma';
colormap='viridis';
colormap='gray';

% 'Clims': Defines the data limits for the color, usefull when one wants to
% see some specific range. The default computes the 5% and 95% percentiles
% of the data and uses that as limit.

clims=[0 0.03];

% 'Savegif': allows to save the plotted figure as an animated gif,
% specified by the given filename.

giffilename='demo5image.gif';

% 'Slice': allows to plot a single slice .Will overwrite the behaviour
% of 'Step'
slice=64;

% Lets go for it

plotImg(imgFDK,'Dim',dimension,'Step',step,'CLims',clims,'Colormap',colormap,'Savegif',giffilename);

% Remember: You can always plot more than 1 image together!
plotImg([thorax imgFDK imgOSSART],'Dim','z')
% Or even the errors!
plotImg([abs(thorax-imgFDK) abs(thorax-imgOSSART)],'Dim',3);