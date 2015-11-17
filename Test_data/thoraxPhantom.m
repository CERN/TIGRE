function img=thoraxPhantom(varargin)
%THORAXPHANTOM returns the thorax phantom 
%
%  IMG=THORAXPHANTOM() returns 128^3 image IMG
%
%  IMG=THORAXPHANTOM(SZ) returns a SZ^3 image IMG if SZ is scalar, or a [SZ(1)
%   SZ(2) SZ(3] image, of SZ is a vector. 


% Deal with inputs
if nargin==0
    sz=[128,128,128];
end
if nargin==1
    if ~isnumeric(varargin{1})
        error('CBCT:thoraxphantom:invalid input','input is not numerical');     
    end
    nsz=max(size(varargin{1}));
    if nsz==2 || nsz>3
        error('CBCT:thoraxphantom:invalid input','input is not 1x1 or 1x3');
    end
    if nsz==1
        sz=[varargin{1},varargin{1},varargin{1}];
    else
        sz=varargin{1};
    end
end
% load data
data=load('img128.mat');
img=double(data.img);

% interpolate data to get desired size
[y, x, z]=...
   ndgrid(linspace(1,size(img,1),sz(1)),...
          linspace(1,size(img,2),sz(2)),...
          linspace(1,size(img,3),sz(3)));
      
imOut=interp3(img,x,y,z,'nearest');
% out!
img=imOut;

end