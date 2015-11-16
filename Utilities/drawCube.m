function drawCube( origin, size,color,alpha,a)
% From
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/235581

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