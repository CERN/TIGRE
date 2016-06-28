function h=plotgeometry(geo,angle)
%PLOTGEOMETRY(GEO,ANGLE) plots a simplified version of the CBCT geometry with the
% given geomerty GEO and angle ANGLE. If angle is nnot Give, 0 will be chosen.
% 
% h=PLOTGEOMETRY(...) will return the figure handle 
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
if nargin<2
    angle=0;
end
angle=angle*180/pi;
%% Figure stuff
h=figure('Name','Cone Beam Compute Tomography geometry');
hold on
title('Current CBCT geometry, in scale')
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gcf, 'color', [1 1 1])

%% CUBE/Image

drawCube(geo.offOrigin,geo.sVoxel,'k',0.05,0);

%% Detector
drawCube([-geo.DSD+geo.DSO; geo.offDetector],[1; geo.sDetector],'r',1,angle)
plotCircle3D([0 0 0],-geo.DSD+geo.DSO);
%% source
p=plot3(geo.DSO,0,0,'.','MarkerSize',30);
rotate(p,[0  0  1],angle,[0 0 0]);
plotCircle3D([0 0 0],geo.DSO);

%% Arrows.
arrow=geo.sVoxel;
%XYZ arrows
try
arrow3d([0 arrow(1)],[0 0],[0 0],.90,5,15,'r');
arrow3d([0 0],[0 arrow(2)],[0 0],.90,5,15,'b');
arrow3d([0 0],[0 0],[0 arrow(3)],.90,5,15,'g');
catch e
    error('CBCT:plotgeometry:Arrow','arrow3D not found, make sure its added to the path');
end
%UV arrows



%%
axis equal;
view(128,26)
% grid on
% axis off;



end
function drawCube( origin, size,color,alpha,a)
% From
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/235581

x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1);
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2);
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3);
for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),color);
    set(h,'facealpha',alpha)
    rotate(h,[0 0 1],a,[0 0 0]);
end

end

function plotCircle3D(center,radius)

theta=0:0.1:2*pi;
v=null([0 0 1]);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'k-.');

end