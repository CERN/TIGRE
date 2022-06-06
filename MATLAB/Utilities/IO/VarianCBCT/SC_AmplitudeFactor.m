function cfactor = SC_AmplitudeFactor(blk, page, edgewt, sccalib)
% Amplitude Factor: ce_i with Edge Response Function (edgewt)
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Author: Yi Du (yi.du@hotmail.com)
% Date: 2021-05-24

%% group number
ngroup = length(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel);

% Amplitude factor groups
cfactor = [];

for ii=1:ngroup
    tmp = sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.ObjectScatterFit;
    % Amplitude Factor
    % unit: mm - > cm
    A = str2double(tmp.A.Text) / 10;
    % unitless
    alpha = str2double(tmp.alpha.Text);
    beta = str2double(tmp.beta.Text);

    % fill holes
    term = (page + eps)./(blk + eps);
    logterm = -log(term);
    logterm(logterm<0) = NaN;
    logterm = single(inpaint_nans(double(logterm), 2));
    
    % amplitude factor wi edge response function as well
    cfactor(:,:,ii) = A .*edgewt .* (term).^(alpha) .* ( logterm ).^(beta);    
end

end

