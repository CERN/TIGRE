function [proj,geo]=flatten_detector(proj,geo)
% This function creates projections on a flat panel detector from a curved
% detector. 2D only (are there any 3D curved detectors?). 

% lets assume we know the arc-length (this has to be input or computed
% somewhow). We assume this is to the end of the pixel, not to the center
% of the last pixel
arclength=45*pi/180;
pixel_arclength=arclength/size(proj,2);


%% useful tools

distance_arc_to_plane=@(alpha)((sqrt(1/(1-sin(alpha).^2))-1)*geo.DSD);

%% Geometric interesting values. 

% dDetector. The smallest projected value will be the central one. 
if mod(size(proj,2),2)
   h=geo.DSD+distance_arc_to_plane(pixel_arclength);
   geo.dDetector=[sqrt(h^2-geo.DSD^2);1];
else
   h=geo.DSD+distance_arc_to_plane(pixel_arclength/2);
   geo.dDetector=[sqrt(h^2-geo.DSD^2)*2;1];
end
% sDetector.
geo.sDetector=[(geo.DSD+distance_arc_to_plane(arclength/2))*2;1];
% nDetector
geo.nDetector=geo.sDetector./geo.dDetector;

%% Interpolation

% We need to define the point in the flat detector
% then move them to the arc, and compute their angular-index. 
% that index will be used to interpolate. 

% step 1: define points
detector=(-geo.sDetector(1)/2+geo.dDetector(1)/2):geo.dDetector(1):(geo.sDetector(1)/2-geo.dDetector(1)/2);
angle=asin(detector.^2./(detector.^2+geo.DSD.^2));
normalized_angle=angle/pixel_arclength;

for ii=1:size(proj,3)
    proj(:,:,ii)=interp1(proj(:,:,ii),normalized_angle);
end

end