function [rotDetector,angles]=coord2angles(geo,det_origin,up_vec,norm_vec)
%COORD2ANGLES converts real world coordinates of the detector to TIGRE angles
%
%
% The detector origin is assumed to have no detector offset!
%
% There are 2 relevant rotations here. The ones around the detector center
% and the ones around the object center. From the origin and 2 vectors we
% can obtain both.

% Input checking:

assert(((geo.DSD-geo.DSO)-sqrt(sum(det_origin(:).^2))<1e-4),'Detector is not at geo.DOD distance from origin.');

assert(dot(up_vec,norm_vec)<1e-4,'Up vector and normal vector are not orthonormal');

% Lets normalize the vectors
up_vec=up_vec/sqrt(sum(up_vec(:).^2));
norm_vec=norm_vec/sqrt(sum(norm_vec(:).^2));

%% first lets obtain local rotation (around detector origin)
% Lets define the vectors we want

non_rotated_norm_vec=-det_origin/sqrt(sum(det_origin(:).^2));

% The up vector is harder, as there are infinte posibilities. We chose the
% one that is angularly closer to the real one. This is done by projecting
% the vector onto the plane
% https://math.stackexchange.com/questions/633181/formula-to-project-a-vector-onto-a-plane

non_rotated_up_vec=up_vec-(dot(up_vec,non_rotated_norm_vec)*non_rotated_norm_vec);

% solve! (Wahba's problem using SVD)

B=up_vec*non_rotated_up_vec.'+norm_vec*non_rotated_norm_vec.'; % this should be 3x3
assert(isequal(size(B),[3 3]),'input vectors are not 3x1, (they are likely 1x3)')

[U,~,V]=svd(B);
M=eye(3); M(3,3)=det(U)*det(V);
R=U*M*V';

% This rotation is roll-pitch-yaw.
% Computing Euler angles from a rotation matrix -Greg Slabaugh

if (abs(R(3,1))-1) < 1e-10
    pitch1=-asin(R(3,1));
    % pitch2=pi-pitch1;
    yaw1=atan2(R(3,2)/cos(pitch1),R(3,3)/cos(pitch1));
    % yaw2=atan2(R(3,2)/cos(pitch2),R(3,3)/cos(pitch2));
    roll1=atan2(R(2,1)/cos(pitch1),R(1,1)/cos(pitch1));
    % roll2=atan2(R(2,1)/cos(pitch2),R(1,1)/cos(pitch2));
else
    roll1=0;
    if sign(R(3,1))==1
        pitch1=-pi/2;
        yaw1=-roll1+atan2(-R(1,2),-R(1,3));
    else
        pitch1=pi/2;
        yaw1=roll1+atan2(R(1,2),R(1,3));
    end
end


rotDetector=[roll1;pitch1;yaw1];

%% Now lets get global rotation around the origin. We need new vectors


origin_norm_vec=det_origin+non_rotated_norm_vec;
origin_up_vec=det_origin+non_rotated_up_vec;

% now we can reuse the variables
non_rotated_norm_vec=[geo.DSO-1;0;0];
non_rotated_up_vec=[geo.DSO;0;1];

B=origin_up_vec*non_rotated_up_vec.'+origin_norm_vec*non_rotated_norm_vec.'; % this should be 3x3
[U,~,V]=svd(B);
M=eye(3); M(3,3)=det(U)*det(V);
R=U*M*V';


if (abs(R(3,3))-1) < 1e-10
    
    angle21=acos(R(3,3));
    angle22=angle21+pi;
    angle31=atan2(R(3,2)/sin(angle21),-R(3,1)/sin(angle21));
    angle32=atan2(R(3,2)/sin(angle22),-R(3,1)/sin(angle22));
    angle11=atan2(R(2,3)/sin(angle21),R(1,3)/sin(angle21));
    angle12=atan2(R(2,3)/sin(angle22),R(1,3)/sin(angle22));
else
    angle31=0;
    if sign(R(3,3))==1 
        angle21=pi/2;
        angle11=-angle31+atan2(R(1,1),R(2,1));
    else
        angle21=-pi/2;
        angle21=pi/2;
        angle11=angle31+atan2(-R(1,1),R(2,1));
    end
end
angles=[angle11;angle21;angle31];