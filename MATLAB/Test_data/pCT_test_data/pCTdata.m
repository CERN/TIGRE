function data=pCTdata()

%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%                     and
%                     https://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=35548
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Stefanie Kaser 
%--------------------------------------------------------------------------

% Load pCT test data (method as in ../MRheadbrain/headPhantom.m by 
% Kyungsang Kim & Ander Biguri )
curr_path=mfilename('fullpath');
data_path=curr_path(1:end-length('/MATLAB/Test_data/pCT_test_data/pCTdata'));
data_path=[data_path '/Common/data/'];
data=load([data_path 'pCT_data.mat']);

end
