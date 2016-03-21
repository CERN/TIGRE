function [p,ellipse]= sheppLogan3D( size,type )
%SHEPPLOGAN3D(SIZE,TYPE) returns the shepp logan phantom, defined by size 
% SIZEand type :
% 
%      'Shepp-Logan'            A test image used widely by researchers in
%                               tomography
%      'Modified Shepp-Logan'   (default) A variant of the Shepp-Logan phantom
%                               in which the contrast is improved for better  
%                               visual perception.
%      'yu-ye-wang'             Another version of the modified Shepp-Logan
%                               phantom from "Katsevich-Type Algorithms for
%                               Variable Radius Spiral Cone-BeamCT"
% 
% Default values are 128^3 and 'Modified Shepp-Logan'

if nargin==0
    sz=[128,128,128];
    type='modified shepp-logan';
else if nargin<2
    type='modified shepp-logan';
end
    
[p,ellipse]=phantom3dAniso(size,type);


end

