function [p,ellipse]= sheppLogan3D( sz,type )
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
if nargin==0
    sz=[128,128,128];
    type='modified shepp-logan';
else
    if nargin<2
    type='modified shepp-logan';
    end
end

warning('This file is NOT under BSD license!')
[p,ellipse]=phantom3dAniso(sz,type);

p=single(p);
% p=[];
end

