function [qualMeas]=Measure_Quality(res_prev,res,QualMeasOpts)
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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Manasavee Lohvithee
%--------------------------------------------------------------------------

%for loop over parameters to get the name 
for ii=1:length(QualMeasOpts)
    opt=QualMeasOpts{ii};
   
    switch opt
        case 'RMSE'
         q=RMSE(res_prev,res);
         
        case 'CC'
         q=CC(res_prev,res); 
         
        case 'MSSIM'
         q=MSSIM(res_prev,res);
         
        case 'UQI'
         q=UQI(res_prev,res);
        case 'error_norm'
         q=im3Dnorm(res_prev-res,'L2');
        
    end
    
    
    qualMeas(ii)=q; 
end



end