%% Demo 10: How to use Image quality measures
%
% This demo demonstrate how to compute all the image quality measures
% by calling the "Measure_Quality.m" function with detailed description.
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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Manasavee Lohvithee
%--------------------------------------------------------------------------
%% Initialize

clear;
close all;
%% Define Geometry
% 
geo=defaultGeometry('nVoxel',[128;128;128]);                     

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi,100);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);

%% "Measure_Quality.m" function
% There are three inputs to the function
% 
% "res_prev"  is a set of real phantom image or reconstructed image from
%             the previous iteration (in case you want to compare the image
%             quality measure after each iteration)
% 
% "res"       is a set of final reconstructed image or reconstructed image
%             from the current iteration (in case you want to compare the image
%             quality measure after each iteration)
%
% "QualMeasOpts"  contains a cell array of desired quality measurement names. 
%                 Example: {'CC','RMSE','MSSIM','UQI'}. At the moment,
%                 there are 4 image quality measures as following:
%
%           'RMSE' : The square root of the mean of the squared differences
%                    of pixel intensities of two images.
%
%                    Taken from "A rigid motion correction method for helical
%                    computed tomography (CT)"
%                    (doi:10.1088/0031-9155/60/5/2047)
%
%           'CC'   : The Pearson correlation coefficient which measures the
%                    linear dependence between two images.
%
%                    Taken from "A rigid motion correction method for helical
%                    computed tomography (CT)"
%                    (doi:10.1088/0031-9155/60/5/2047)
%
%           'MSSIM': The mean structural similarity index, a measure of the
%                    similarity of two images in terms of luminance,
%                    contrast and structure that is designed to provide a
%                    good approximation of perceptual image quality.
%
%                    Taken from "A rigid motion correction method for helical
%                    computed tomography (CT)"
%                    (doi:10.1088/0031-9155/60/5/2047)
%
%           'UQI' : The universal quality index to evaluate the degree of
%                   similarity between the reconstructed and phantom images
%                   for chosen ROIs. Its value ranges from zero to one.
%                   A UQI value closer to one suggests better similarity to true image.
%
%                   Taken from "Few-view cone-beam CT reconstruction with
%                   deformed prior image" (doi: 10.1118/1.4901265)
%
%
%
%

%% Define a cell array of image quality measures 
qualmeas={'RMSE','CC','MSSIM','UQI'};

blcks=22;

%% Reconstruct image using SIRT, SART and OS-SART

[imgSIRT,errL2SIRT,qualitySIRT]=SIRT(noise_projections,geo,angles,30,...
                            'QualMeas',qualmeas);
[imgSART,errL2SART,qualitySART]=SART(noise_projections,geo,angles,30,...
                            'QualMeas',qualmeas);
[imgOSSART,errL2OSSART,qualityOSSART]=OS_SART(noise_projections,geo,angles,30,...
                            'QualMeas',qualmeas,'BlockSize',blcks);
                             
                        
                        
%% Plot the changes of image quality measures over iterations

%%RMSE plot
figure
plot([qualitySIRT(1,:);[qualitySART(1,:) nan(1,length(qualityOSSART)-length(qualitySART))];qualityOSSART(1,:)]');
title('Evolution of RMSE per iteration')
legend('SIRT','SART','OS-SART')
                   
%%CC plot
figure
plot([qualitySIRT(2,:);[qualitySART(2,:) nan(1,length(qualityOSSART)-length(qualitySART))];qualityOSSART(2,:)]');
title('Evolution of correlation coefficient (CC) per iteration')
legend('SIRT','SART','OS-SART')

%%MSSIM plot
figure
plot([qualitySIRT(3,:);[qualitySART(3,:) nan(1,length(qualityOSSART)-length(qualitySART))];qualityOSSART(3,:)]');
title('Evolution of mean structural similarity index per iteration')
legend('SIRT','SART','OS-SART')
                        
%%UQI plot
figure
plot([qualitySIRT(4,:);[qualitySART(4,:) nan(1,length(qualityOSSART)-length(qualitySART))];qualityOSSART(4,:)]');
title('Evolution of universal quality index per iteration')
legend('SIRT','SART','OS-SART')

%% Compute the parameters for the target image
% Knowing change of the parameters per iteration can be interesting, but
% the main use of them would be to compare a reconstructed image to a given
% known target image. 

disp('UQI for SART,   SIRT,   OS-SART')
disp(['       ',num2str(UQI(head,imgSART)),' ', num2str(UQI(head,imgSIRT)),' ' ,num2str(UQI(head,imgOSSART))])
disp('RMSE for SART,   SIRT,   OS-SART')
disp(['       ',num2str(RMSE(head,imgSART)),' ', num2str(RMSE(head,imgSIRT)),' ' ,num2str(RMSE(head,imgOSSART))])