%Creating CUDA accelerated test projection

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
% Coded by:           Stefanie Kaser 
%--------------------------------------------------------------------------

%Geometry initialization
% pCT geometry definition
E_INIT = single(100);
PIXEL_SIZE = single(1);
POS_DET_IN = single(-100);
POS_DET_OUT = single(100);
DETECTOR_SIZE_X = single(500);
DETECTOR_SIZE_Y = single(500);
HULL_PARAMS = single([0; 0; 0; 0]); 

numOfProtons = 10000000;
%------------------ Data aquisition --------------------------

tic

posIn = single(zeros(2*numOfProtons, 1));
posOut = single(zeros(2*numOfProtons, 1));
dirIn = single(zeros(2*numOfProtons, 1));
dirOut = single(zeros(2*numOfProtons, 1));
Wepl = single(zeros(numOfProtons, 1)); 

for k = 1:numOfProtons
    posIn(k) = single(-(DETECTOR_SIZE_X/2)+(DETECTOR_SIZE_X)*rand);
    posIn(k+numOfProtons) = single(-(DETECTOR_SIZE_Y/2)+(DETECTOR_SIZE_Y)*rand);
    posOut(k) = single(posIn(k));
    posOut(k+numOfProtons) = single(posIn(k+numOfProtons));
    if sqrt(posOut(k)^2 + posOut(k+numOfProtons)^2) < DETECTOR_SIZE_X/4
        Wepl(k) = single(0);
    else
        Wepl(k) = single(100);
    end
    
end
toc

tic
proj = pCTCubicSpline_mex(posIn, posOut, dirIn, dirOut, Wepl, PIXEL_SIZE, DETECTOR_SIZE_X, DETECTOR_SIZE_Y, POS_DET_IN, POS_DET_OUT, E_INIT, HULL_PARAMS);
toc
imshow(proj);