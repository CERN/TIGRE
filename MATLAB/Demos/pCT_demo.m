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
geo.dDetector = [2; 2];
geo.DSID = 500;
geo.DSO = 600;
geo.DSD = 700;
geo.hull = [0; 0; 0; 0];
geo.sDetector = [500; 500];
geo.mode = 'parallel';
% pCT energy definition
eIn = single(100);
numOfProtons = 1000000;
%------------------ Data aquisition --------------------------

posIn = single(zeros(2*numOfProtons, 1));
posOut = single(zeros(2*numOfProtons, 1));
dirIn = single(zeros(2*numOfProtons, 1));
dirOut = single(zeros(2*numOfProtons, 1));
Wepl = single(zeros(numOfProtons, 1)); 


for k = 1:numOfProtons
    posIn(k) = single((-1*geo.sDetector(1)/2)+(geo.sDetector(1)*rand));
    posIn(k+numOfProtons) = single((-1*geo.sDetector(2)/2)+(geo.sDetector(2)*rand));
    posOut(k) = single(posIn(k));
    posOut(k+numOfProtons) = single(posIn(k+numOfProtons));
    if sqrt(posOut(k)^2 + posOut(k+numOfProtons)^2) < geo.sDetector(1)/4
        Wepl(k) = single(100);
    else
        Wepl(k) = single(0);
    end
    
end

tic
proj = pCTCubicSpline_mex(posIn, posOut, dirIn, dirOut, Wepl, eIn, geo);
toc
imshow(proj);