%% Defining Geometry
% Distances
geo.DSD = 750;                             % Distance Source Detector      (mm)
geo.DSO = 380;                             % Distance Source Origin        (mm)
% Detector parameters
geo.nDetector=[650; 650];					% number of pixels              (px)
geo.dDetector=[0.24; 0.24]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[256;256;256];                   % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
                                            % These two can be also defined

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi-2*pi/50,50);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');

niter=50;

%% Algorithm
% 1. FDK with hann filter.
imgFDK_hann=FDK(projections,geo,angles,'filter','hann');

% 2. FDK with hann filter.
imgFDK_ramlak=FDK(projections,geo,angles,'filter','ram-lak');

% 3. SIRT
lambda=1;
lambdared=0.999;
initmode='none';
verbose=true;
qualmeas={'RMSE'};
[imgSIRT,errL2SIRT,qualitySIRT]=SIRT(projections,geo,angles,20,...
                            'lambda',lambda,'lambda_red',lambdared,'verbose',verbose,'QualMeas',qualmeas);

% 4. SART                        
[imgSART,errL2SART,qualitySART]=SART(projections,geo,angles,20,...
                            'lambda',lambda,'lambda_red',lambdared,'verbose',verbose,'QualMeas',qualmeas);

%5. OS-SART
blcks=8;
order='angularDistance';
[imgOSSART,errL2OSSART,qualityOSSART]=OS_SART(projections,geo,angles,20,...
                            'lambda',lambda,'lambda_red',lambdared,'verbose',verbose,'QualMeas',qualmeas,...
                             'BlockSize',blcks,'OrderStrategy',order);
                         
%6. ASD-POCS
% [imgOSSART_ASDPOCS,errL2OSSART_ASDPOCS]=OS_SART(projections,geo,angles,60);
epsilon=errL2OSSART(end);
alpha=0.002;
ng=25;
lambda=1;
lambdared=0.9999;
alpha_red=0.95;
ratio=0.94;
verb=true;

imgASDPOCS=ASD_POCS(projections,geo,angles,50,...
                    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb); % less important.

%7. OS-ASD-POCS
imgOSASDPOCS=OS_ASD_POCS(projections,geo,angles,50,...
                    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb,...% less important.
                      'BlockSize',size(angles,2)/10,'OrderStrategy','angularDistance'); %OSC options

% Plotting Images
plotImg([imgFDK_hann, imgFDK_ramlak, imgSIRT, imgSART, imgOSSART, imgASDPOCS, imgOSASDPOCS],'Dim','Z', 'Slice',92:110,'clims',[0 2]);