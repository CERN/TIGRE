function drawCube( origin, size,color,alpha,a)
% From
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/235581
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
x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1);
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2);
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3);
arad=a*pi/180;
R=[cos(arad) -sin(arad) 0; sin(arad) cos(arad) 0; 0 0 1];

for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),color);
    set(h,'facealpha',alpha)
    rotate(h,[0 0 1],a,[0 0 0]);
end

end