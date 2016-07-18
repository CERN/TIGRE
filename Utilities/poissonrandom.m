function km=poissonrandom(lambda)

% Donald Knuth's poisson random generator

% This sucks bananas because its O(lambda) and we want lambda~=10000...
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

km=zeros(size(lambda));
for ii=1:numel(lambda)
    
    L=exp(-lambda(ii));
    k=1;
    p=1*rand;
    
    while p>L
        k=k+1;
        p=p*rand;
    end
    km(ii)=k-1;
end
end