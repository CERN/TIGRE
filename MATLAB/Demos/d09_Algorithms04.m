%% DEMO 9:  Algorithms 04. Total Variation minimization algorithms
%
%
%  This demo presents the Total variation algorithms in TIGRE. Total
%  variation algorithms try to minimize the variation (gradient) of the
%  image, assuming its piecewise smooth, as most things in nature are (i.e.
%  human body). 
%
% This set of algorithms is specially good performing when the noise is
% very big or the number of projections is small, however, they require more
% computational time and memory than the other algorithms to run.
% 
%
% This Demo is not fully complete, as several algorithms are not shown
% here, yet are available in TIGRE, namely:
%
% - AwASD_POCS:  Edge preserving ASD_POCS (Adaptative weigthed).
% - OS_AwASD_POCS: OS-version of the previous algorithm
% - PCSD: A version of ASD_POCS that heuristically select some of the
%         parameters, particularly epsilon (maxL2norm)
% - AwPCSD: Edge preserving version of the previous algorithm
% - OS_AwPCSD: block-wise version of the previous
%
% You can use these algorithms in a very similar manner as below examples.
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
geo=defaultGeometry('nVoxel',[128;128;128]);                     

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi-2*pi/30,30);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);

%% Lets create a OS-SART test for comparison
[imgOSSART,errL2OSSART]=OS_SART(noise_projections,geo,angles,60);

%% Total Variation algorithms
%
%  ASD-POCS: Adaptative Steeppest Descent-Projection On Convex Subsets
% Often called POCS-TV
%==========================================================================
%==========================================================================
%  ASD-POCS minimizes At-B and the TV norm separately in each iteration,
%  i.e. reconstructs the image first and then reduces the TV norm, every
%  iteration. As the other algorithms the mandatory inputs are projections,
%  geometry, angles and maximum iterations.
%
% ASD-POCS has a veriety of optional arguments, and some of them are crucial
% to determine the behaviour of the algorithm. The advantage of ASD-POCS is
% the power to create good images from bad data, but it needs a lot of
% tunning. 
%
%  Optional parameters that are very relevant:
% ----------------------------------------------
%    'maxL2err'    Maximum L2 error to accept an image as valid. This
%                  parameter is crucial for the algorithm, determines at
%                  what point an image should not be updated further.
%                  Default is 20% of the FDK L2 norm.
%                  
% its called epsilon in the paper
epsilon=im3Dnorm(Ax(FDK(noise_projections,geo,angles),geo,angles)-noise_projections,'L2')*0.15;
%   'alpha':       Defines the TV hyperparameter. default is 0.002. 
%                  However the paper mentions 0.2 as good choice
alpha=0.002;

%   'TViter':      Defines the amount of TV iterations performed per SART
%                  iteration. Default is 20

ng=25;

% Other optional parameters
% ----------------------------------------------
%   'lambda':      Sets the value of the hyperparameter for the SART iterations. 
%                  Default is 1
%
%   'lambdared':   Reduction of lambda Every iteration
%                  lambda=lambdared*lambda. Default is 0.99
%
lambda=1;
lambdared=0.9999;


%   'alpha_red':   Defines the reduction rate of the TV hyperparameter
alpha_red=0.95;

%   'Ratio':       The maximum allowed image/TV update ration. If the TV 
%                  update changes the image more than this, the parameter
%                  will be reduced.default is 0.95
ratio=0.94;

%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.

verb=true;

imgASDPOCS=ASD_POCS(noise_projections,geo,angles,50,...
                    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb); % less important.



                 
%  OS_ASD_POCS: Odered Subset-TV algorithm
%==========================================================================
%==========================================================================
%
% The logical next step to imporce ASD-POCS is substituting SART with a
% faster algorithm, such as OS-SART
%
% The parameters are the same as in ASD-POCS, but also have 'BlockSize' and
% @OrderStrategy', taken from OS-SART

imgOSASDPOCS=OS_ASD_POCS(noise_projections,geo,angles,50,...
                    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb,...% less important.
                      'BlockSize',size(angles,2)/10,'OrderStrategy','angularDistance'); %OSC options
           
                
           
%  B-ASD-POCS-beta 
%==========================================================================
%==========================================================================
% Is another version of ASD-POCS that has uses Bregman iteration, i.e.
% updates the data vector every some ASD-POCS iteration, bringing closer
% the data to the actual image. According to the original article, this
% accelerates convergence rates, giving a better image in the same amount
% of iteratiosn.
%
% It has all the inputs from ASD-POCS, and has 3 extra:
%   'beta'         hyperparameter controling the Bragman update. default=1
% 
%   'beta_red'     reduction of the beta hyperparameter. default =0.75
% 
%   'bregman_iter' amount of global bregman iterations. This will define
%                  how often the bregman iteration is executed. It has to
%                  be smaller than the number of iterations.

%Note that the number of iteration for TV has changed

imgBASDPOCSbeta=B_ASD_POCS_beta(noise_projections,geo,angles,50,...
                    'TViter',40,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb,... % less important.
                      'beta',0.5,'beta_red',0.7,'bregman_iter',10); % bregman options
                  
%   SART-TV 
%==========================================================================
%==========================================================================      
%
%   This implementation differs more from the others, as it minimizes the
%   ROF model, i.e. when minimizing the total variation of the image, it
%   also takes into account the information of the image. If only the total
%   variation minimization step was run in the rest of the algorithms, the
%   result would be a flat image (as that is the minimum total variation
%   image), altertatively, the ROF model enforces the image not to change too
%   much.
%   
%   This algirths performs better with more projections, but noisy data, by
%   enforncing the TV very little
%  
%   The optional parameters are for the total variatiot part of the image:
%
%
%  
%   'TViter'       amoutn of iteration in theTV step. Default 50
% 
%   'TVlambda'     hyperparameter in TV iteration. IT gives the ratio of
%                  importance of the image vs the minimum total variation.
%                  default is 15. Lower means more TV denoising.
% 
                  
imgSARTTV=SART_TV(noise_projections,geo,angles,50,'TViter',100,'TVlambda',50);           

 %% Lets visualize the results
% Notice the smoother images due to TV regularization.
%
%     thorax              OS-SART           ASD-POCS         
%    
%     OSC-TV             B-ASD-POCS-beta   SART-TV

plotImg([ imgOSASDPOCS imgBASDPOCSbeta imgSARTTV; head imgOSSART  imgASDPOCS ] ,'Dim','Z','Step',2,'clims',[0 1])
 % error

plotImg(abs([ head-imgOSASDPOCS head-imgBASDPOCSbeta head-imgSARTTV;head-head head-imgOSSART  head-imgASDPOCS ]) ,'Dim','Z','Slice',64)




%%

% Obligatory XKCD reference: https://xkcd.com/833/



