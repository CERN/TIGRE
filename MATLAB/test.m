%% Init 
clear;
close all;
%% Test parameters
test_case = 3

if test_case == 0
    % Both even
    nVoxelZ=10;
    nDetectorV=10;
elseif test_case == 1
    % nVoxel Even, nDetec Odd
    nVoxelZ=10;
    nDetectorV=11;
elseif test_case == 2
    % nVoxel Odd, nDetec Even
    nVoxelZ=11;
    nDetectorV=10;
else
    % nVoxel Odd, nDetec Odd
    nVoxelZ=11;
    nDetectorV=11;
end

nVoxelX = nVoxelZ;
nVoxelY = nVoxelZ;
nDetectorU = nDetectorV;

nangles=10;

mode='parallel';  % 'parallel'
ptype= 'Siddon'; % 'interpolated' % 'ray-voxel'
%%
geo=defaultGeometry('mode',mode);
geo.DSD = 2000;
geo.DSO = 1000;
geo.nVoxel = [nVoxelX;nVoxelY;nVoxelZ];
geo.dVoxel = [1.0; 1.0; 1.0];
geo.sVoxel = geo.nVoxel .* geo.dVoxel;
geo.dDetector = [1;1] ;              % size of each pixel            (mm)
geo.nDetector = [nDetectorU; nDetectorV];     % has to be true 3D for python version
geo.sDetector = geo.dDetector .* geo.nDetector;
geo.offOrigin = [0;0;0];
geo.offDetector = [0.0; 0.0] ;
geo.rotDetector = [0; 0; 0];
geo.accuracy=0.01;

%%
angles = linspace(0, 2 * pi - 2*pi/nangles, nangles);

wire = zeros(geo.nVoxel', 'single');
wire(int32((geo.nVoxel(1))/2):int32((geo.nVoxel(1)+1)/2), int32((geo.nVoxel(2))/2):int32((geo.nVoxel(2)+1)/2), :) = 1;
%%
proj = Ax(wire,geo,angles, ptype);
%% Plot
% cellLeg = {};
% for iv =1:nDetectorV
%     plot(squeeze(proj(iv, int32(nDetectorU/2), :)));
%     cellLeg{iv} = sprintf('V=%d', iv);
%     hold on
% end
% legend(cellLeg);
% hold off
plotProj(proj, angles)
