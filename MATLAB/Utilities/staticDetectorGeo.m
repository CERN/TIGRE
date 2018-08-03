function geo=staticDetectorGeo(geo,angles)
% This fucntion computes the translation and rotations needed in the system
% to desribe a detector that does not move when the source moves.
R=(geo.DSD-geo.DSO); %Radious of rotation of detector
geo.DSD=geo.DSD-(R-R*cos(angles));
geo.offDetector=[R*sin(angles); zeros(1,size(angles,2))];
geo.rotDetector=[zeros(1,size(angles,2));zeros(1,size(angles,2));-angles];


%% Hum, I think this is right, but I have been told that the otherone is better.
% % Only works on (-90,90 deg)
% R=(geo.DSD-geo.DSO); 
% geo.DSD=geo.DSD+(R./sin(pi/2-angles)-R);
% geo.offDetector=[R*sin(angles)./sin(pi/2+angles); zeros(1,size(angles,2))];
% geo.rotDetector=[zeros(1,size(angles,2));zeros(1,size(angles,2));-angles];
