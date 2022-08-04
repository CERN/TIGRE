%% Demo 21: Loading scanner data to TIGRE. 
%
% This demo will demostrate the options for loading scanner data into
% TIGRE.
%
%   Supported Manufacturers:
%
%       Varian
%
%       Nikon
%
%       YXLON     
%
%       Xradia (Zeiss)
%
%       Philips (in theory, any DICOM, but only tested in Philips Allura)
%
% Currently we have instructions for generic, Nikon (micro-CT) and Varian
% scanners. 
%
% We are always looking to expand this, if you have code for a scanner that
% is not supported by TIGRE and are allowed to share it, please do, we'd
% like to have as many as possible. 

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
%%

%% Varan onboard CBCT.

% If you have a Varian scanner that saves data in XIM format, the following
% code will be able to load your dataset. 
%
% If you have a recosntructed image, it will take the geometry from there,
% but if you don't, it will estimate appropiate geometry. You can always
% change the image parameters, but the detector parameters should stay the
% same unless you change the projections. 

datafolder='~/your_data_path/varian/2020-01-01_123456/';
[proj,geo, angles ] = VarianDataLoader(datafolder);
[proj,geo, angles ] = VarianDataLoader(datafolder, false); %remove motion lag correction.

% You can directly call reconstruction code now:

img=OS_SART(proj,geo,angles,100);
img=FDK(proj,geo,angles);

%% Nikon micro-CT

% Similarly, a Nikon dataset can be loaded with the following code:

datafolder='~/your_data_path/Nikon/Sample_name/';
[proj,geo, angles ] = NikonDataLoader(datafolder);

% as micro-CT datasets are large, optional arguments for loading partial
% amount of data is available:

% load equidistant angles, but only few:
[proj,geo, angles ] = NikonDataLoader(datafolder,'sampling','equidistant','num_angles',150);
% load every X angles (10)
[proj,geo, angles ] = NikonDataLoader(datafolder,'sampling','step','sampling_step',10);
% load first X angles (1000)
[proj,geo, angles ] = NikonDataLoader(datafolder,'sampling','continuous','num_angles',1000);

% You can directly call reconstruction code now:

img=OS_SART(proj,geo,angles,100);
img=FDK(proj,geo,angles);


%% YXLON

% You can replace NikonDataLoader in the avobe code by YXLONDataLoader().
% the rest of the functionality is the same. 

datafolder='~/your_data_path/YXLON/Sample_name/';
[proj,geo, angles ] = YXLONDataLoader(datafolder);

%% DICOM data (only tested on Philips Allura)

% As with the other loaders, this can be simply done with:

datafolder='~/your_data_path/Dicom/Some_folder/';
[proj,geo, angles,dicomhdr] = dicomDataLoader(datafolder);

%This also returns the headers for DICOM files, so you can inspect them
%yourself. Its been almost untested aside from few datasets, please do
%contact us to help us increase support to more devices/DICOM attributes if
%it does not work for your data. 

% You can directly call reconstruction code now:

img=OS_SART(proj,geo,angles,100);
img=FDK(proj,geo,angles);
%% Generic

% It is possible that your scanner is not currently supported, or that it
% simply does not have any way of storing the information (Zeiss Xradia
% does not store anything but the projecions)

% if this is the case, the general way of tackling the problem is:

% Step 1: define your geometry and angles
    
% Step 2: Load projections:
[num_files, filenames] % Figure out a way of finding this. Often dir([path,'/*.extension'])
% prealocate
proj = single(zeros(geo.nDetector(1),geo.nDetector(2),length(angles)));
% load files
for ii=1:num_files
    proj(:,:,ii)=single(imread(filenames(ii)));
end

% Step 3: validate

plotProj(proj,angles)
% you need to make sure that:
%     1) white=metal/bone/high density and black=air.
%         If this is not the case
proj=-log(proj/(max(proj(:))+1)); % Beer-Lanbert law
%     2) rotation happens left-right instead of top-bottom.
%        If its top bottom:
proj=permute(proj,[2,1,3]);

% Step 4: test
imgfdk=FDK(proj,geo,angles);
plotImg(imgFDK,'dim','z');
% If this does not look good, possible things to check:
% - Are the angles in the right direction? maybe they need to be inverted. 
% - Is your geometry properly defined? mistake on geometry will be
%   detrimental to image quality
% - if microCT: halos around the features? COR correction needed. 
