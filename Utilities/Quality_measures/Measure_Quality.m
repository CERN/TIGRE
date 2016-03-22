function [qualMeas]=Measure_Quality(res_prev,res,QualMeasOpts)


%for loop over parameters to get the name 
for ii=1:length(QualMeasOpts)
    opt=QualMeasOpts{ii};
   
    switch opt
        %%%%%%%%RMSE
        case 'RMSE'
         q=RMSE(res_prev,res);
         
        %%%%%%CC
        case 'CC'
         q=CC(res_prev,res); 
         
        %%%%%%MSSIM
        case 'MSSIM'
         q=MSSIM(res_prev,res);
         
        case 'UQI'
         q=UQI(res_prev,res);
         
        
    end
    
    
    qualMeas(ii)=q; 
end



end