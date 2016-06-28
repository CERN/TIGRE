function [p,ellipse]=phantom3dAniso(varargin)

%PHANTOM3D Three-dimensional analogue of MATLAB Shepp-Logan phantom
%   P = PHANTOM3D(DEF,N) generates a 3D head phantom that can   
%   be used to test 3-D reconstruction algorithms.
%
%   DEF is a string that specifies the type of head phantom to generate.
%   Valid values are: 
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
%   N specifies the 3D grid size of P
%   If N is a scalar, P will have isotropic size [N, N, N]
%   If N is a 3-vector, P will have size [N(1) N(2) N(3)]
%   If you omit the argument, N defaults to [64 64 64].
% 
%   P = PHANTOM3D(E,N) generates a user-defined phantom, where each row
%   of the matrix E specifies an ellipsoid in the image.  E has ten columns,
%   with each column containing a different parameter for the ellipsoids:
%   
%     Column 1:  A      the additive intensity value of the ellipsoid
%     Column 2:  a      the length of the x semi-axis of the ellipsoid 
%     Column 3:  b      the length of the y semi-axis of the ellipsoid
%     Column 4:  c      the length of the z semi-axis of the ellipsoid
%     Column 5:  x0     the x-coordinate of the center of the ellipsoid
%     Column 6:  y0     the y-coordinate of the center of the ellipsoid
%     Column 7:  z0     the z-coordinate of the center of the ellipsoid
%     Column 8:  phi    phi Euler angle (in degrees) (rotation about z-axis)
%     Column 9:  theta  theta Euler angle (in degrees) (rotation about x-axis)
%     Column 10: psi    psi Euler angle (in degrees) (rotation about z-axis)
%
%   For purposes of generating the phantom, the domains for the x-, y-, and 
%   z-axes span [-1,1].  Columns 2 through 7 must be specified in terms
%   of this range.
%
%   [P,E] = PHANTOM3D(...) returns the matrix E used to generate the phantom.
%
%   Class Support
%   -------------
%   All inputs must be of class double.  All outputs are of class double.
%
%   Remarks
%   -------
%   For any given voxel in the output image, the voxel's value is equal to the
%   sum of the additive intensity values of all ellipsoids that the voxel is a 
%   part of.  If a voxel is not part of any ellipsoid, its value is 0.  
%
%   The additive intensity value A for an ellipsoid can be positive or negative;
%   if it is negative, the ellipsoid will be darker than the surrounding pixels.
%   Note that, depending on the values of A, some voxels may have values outside
%   the range [0,1].
%    
%   Example
%   -------
%        ph = phantom3d(128);
%        figure, imshow(squeeze(ph(64,:,:)))
%
%   Copyright 2005 Matthias Christian Schabel (matthias @ stanfordalumni . org)
%   University of Utah Department of Radiology
%   Utah Center for Advanced Imaging Research
%   729 Arapeen Drive
%   Salt Lake City, UT 84108-1218
%
%   This code is released under the Gnu Public License (GPL). For more information, 
%   see : http://www.gnu.org/copyleft/gpl.html
%
%   Portions of this code are based on phantom.m, copyrighted by the Mathworks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modification May 25, 2015, by Patrick J. Bolan, University of Minnesota
% Added support for anisotropic phantom sizes: the phantom size can now be 
% a vector. 
%


[ellipse,n] = parse_inputs(varargin{:});

nx = n(1); ny = n(2); nz = n(3);
p = zeros([nx ny nz]);

rngx =  ( (0:nx-1)-(nx-1)/2 ) / ((nx-1)/2); 
rngy =  ( (0:ny-1)-(ny-1)/2 ) / ((ny-1)/2); 
rngz =  ( (0:nz-1)-(nz-1)/2 ) / ((nz-1)/2); 

% PJB: Note the swap of the x and y with meshgrid parameters. 
%[x,y,z] = meshgrid(rngx,rngy,rngz);
[x,y,z] = meshgrid(rngy,rngx,rngz);

coord = [flatten(x); flatten(y); flatten(z)];

p = flatten(p);

for k = 1:size(ellipse,1)    
   A = ellipse(k,1);            % Amplitude change for this ellipsoid
   asq = ellipse(k,2)^2;        % a^2
   bsq = ellipse(k,3)^2;        % b^2
   csq = ellipse(k,4)^2;        % c^2
   x0 = ellipse(k,5);           % x offset
   y0 = ellipse(k,6);           % y offset
   z0 = ellipse(k,7);           % z offset
   phi = ellipse(k,8)*pi/180;   % first Euler angle in radians
   theta = ellipse(k,9)*pi/180; % second Euler angle in radians
   psi = ellipse(k,10)*pi/180;  % third Euler angle in radians
   
   cphi = cos(phi);
   sphi = sin(phi);
   ctheta = cos(theta);
   stheta = sin(theta);
   cpsi = cos(psi);
   spsi = sin(psi);
   
   % Euler rotation matrix
   alpha = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
            -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
            stheta*sphi                  -stheta*cphi                ctheta];        
   
   % rotated ellipsoid coordinates
   coordp = alpha*coord;
   
   idx = find((coordp(1,:)-x0).^2./asq + (coordp(2,:)-y0).^2./bsq + (coordp(3,:)-z0).^2./csq <= 1);
   p(idx) = p(idx) + A;
end

%p = reshape(p,[nx ny nz]);
p = reshape(p, [nx ny nz]);

return;


function out = flatten(in)

out = reshape(in,[1 numel(in)]);

return;
   
   
function [e,n] = parse_inputs(varargin)
%  e is the m-by-10 array which defines ellipsoids
%  n is a 3-vector with the size of the phantom brain image, [nx ny nz]

n = [64 64 64];     % The default size
e = [];
defaults = {'shepp-logan', 'modified shepp-logan', 'yu-ye-wang'};

for i=1:nargin
   if ischar(varargin{i})         % Look for a default phantom
      def = varargin{i};
      idx = strcmpi(def, defaults);
      if isempty(idx)
         eid = sprintf('Images:%s:unknownPhantom',mfilename);
         msg = 'Unknown default phantom selected.';
         error(eid,'%s',msg);
      end
      switch defaults{idx}
      case 'shepp-logan'
         e = shepp_logan;
      case 'modified shepp-logan'
         e = modified_shepp_logan;
      case 'yu-ye-wang'
         e = yu_ye_wang;
      end
   elseif numel(varargin{i})==1 
      n = [varargin{i} varargin{i} varargin{i}];   % a scalar is the image size
   elseif numel(varargin{i})==3 
      siz = varargin{i};
      n = [siz(1) siz(2) siz(3)]; % 3 integers specify image dimensions   
   elseif ndims(varargin{i})==2 && size(varargin{i},2)==10 
      e = varargin{i};            % user specified phantom
   else
      eid = sprintf('Images:%s:invalidInputArgs',mfilename);
      msg = 'Invalid input arguments.';
      error(eid,'%s',msg);
   end
end

% ellipse is not yet defined
if isempty(e)                    
   e = modified_shepp_logan;
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Default head phantoms:   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function e = shepp_logan

e = modified_shepp_logan;
e(:,1) = [1 -.98 -.02 -.02 .01 .01 .01 .01 .01 .01];

return;

      
function e = modified_shepp_logan
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .810      0       0       0      0      0      0
        -.8  .6624  .874  .780      0  -.0184       0      0      0      0
        -.2  .1100  .310  .220    .22       0       0    -18      0     10
        -.2  .1600  .410  .280   -.22       0       0     18      0     10
         .1  .2100  .250  .410      0     .35    -.15      0      0      0
         .1  .0460  .046  .050      0      .1     .25      0      0      0
         .1  .0460  .046  .050      0     -.1     .25      0      0      0
         .1  .0460  .023  .050   -.08   -.605       0      0      0      0
         .1  .0230  .023  .020      0   -.606       0      0      0      0
         .1  .0230  .046  .020    .06   -.605       0      0      0      0 ];
       
return;
          

function e = yu_ye_wang
%
%   Yu H, Ye Y, Wang G, Katsevich-Type Algorithms for Variable Radius Spiral Cone-Beam CT
%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .900      0       0       0      0      0      0
        -.8  .6624  .874  .880      0       0       0      0      0      0
        -.2  .4100  .160  .210   -.22       0    -.25    108      0      0
        -.2  .3100  .110  .220    .22       0    -.25     72      0      0
         .2  .2100  .250  .500      0     .35    -.25      0      0      0
         .2  .0460  .046  .046      0      .1    -.25      0      0      0
         .1  .0460  .023  .020   -.08    -.65    -.25      0      0      0
         .1  .0460  .023  .020    .06    -.65    -.25     90      0      0
         .2  .0560  .040  .100    .06   -.105    .625     90      0      0
        -.2  .0560  .056  .100      0    .100    .625      0      0      0 ];
       
return;
        
             