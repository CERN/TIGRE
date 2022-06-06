function edgewt = SC_EdgeResponse(thickness)
% Ideal model: k =  EdgeCoef0, (0.4), EdgeRange = EdgeRangeMM, (80)
% i.e., y = k*x + 0.78
% Herein, we used an approximation method to model the edge response
% function
% Date: 2021-05-24


% Emperical Thickness Threshold for Edge Segementation
edgewt = double(imbinarize(thickness, 50));
tmpmask = edgewt;
h = fspecial('average', [25, 25]);
for ii = 1:5
    edgewt = imfilter(edgewt, h);
end        
tmp = tmpmask.*edgewt;
tmp(tmp==0) = NaN;
tmp = rescale(tmp, 0.6, 1);
tmp(isnan(tmp)) = 0;
edgewt = tmp;

end

