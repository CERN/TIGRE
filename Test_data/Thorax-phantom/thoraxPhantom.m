function img=thoraxPhantom(varargin)
%THORAXPHANTOM returns the thorax phantom 
%
%  IMG=THORAXPHANTOM() returns 128^3 image IMG
%
%  IMG=THORAXPHANTOM(SZ) returns a SZ^3 image IMG if SZ is scalar, or a [SZ(1)
%   SZ(2) SZ(3] image, of SZ is a vector. 
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
%                     and
%                     https://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=35548
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Kyungsang Kim, modified by Ander Biguri 
%--------------------------------------------------------------------------

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
img=data.img;

% interpolate data to get desired size
[y, x, z]=...
   ndgrid(linspace(1,size(img,1),sz(1)),...
          linspace(1,size(img,2),sz(2)),...
          linspace(1,size(img,3),sz(3)));
      
imOut=interp3(img,x,y,z,'nearest');
% out!
img=imOut;

end