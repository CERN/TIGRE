function tf = haveInet()
%HAVEINET returns true if the owner has internet or false if it hasent
% Taken from
% https://stackoverflow.com/questions/19557118/internet-connection-status-using-matlab
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
  tf = false;
  try
    java.net.InetAddress.getByName('www.google.com');
    tf = true;
  end
end