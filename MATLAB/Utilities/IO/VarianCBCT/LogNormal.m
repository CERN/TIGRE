function proj = LogNormal(proj, angles, airnorm, Blk, Sec, BlkAirNorm,gpuids)
% Log Normalization: 
% Calculate Logrithmic Projections with AirNorm
% Tested on TrueBeam 2.5 and 2.7
% Date: 2021-03-23
% Author: Yi Du (yi.du@hotmail.com)

%% TB version: CPU
%{
%
tic
% Version = 2.5
if(isempty(Sec))
    for ii = 1:length(angles)
        CF = airnorm(ii)/BlkAirNorm;
        proj(:,:,ii) = log(CF*Blk./(proj(:,:,ii) + eps) + eps);
    end
% Version = 2.7    
else
    % GantryRtn = KvSourceRtn - 90;
    angles = angles - 90;

    % Flip for interpolation
    if(Sec(2) - Sec(1) <0)
        Sec = flip(Sec);
        BlkAirNorm = flip(BlkAirNorm);
        Blk = flip(Blk, 3);
    end    
    
    % interpolation weights
    for ii = 1:length(angles)
        [left, weights] = interp_weight(angles(ii), Sec);
        interp_blk = weights(1) * Blk(:,:,left) + weights(2) * Blk(:,:,left+1);
        % Correction factor
        CF = airnorm(ii)/(0.5*BlkAirNorm(left) + 0.5*BlkAirNorm(left+1));
        proj(:,:,ii) = log(CF*interp_blk./(proj(:,:,ii)+eps) + eps);
    end
end
toc
%}
%% TB version: GPU
gpuDevice(gpuids.devices(0)+1);
gproj = gpuArray(single(proj));
gBlk = gpuArray(Blk);
% Version = 2.5
if(isempty(Sec))
    for ii = 1:length(angles)
        CF = airnorm(ii)/BlkAirNorm;
        gproj(:,:,ii) = log(CF*gBlk./(gproj(:,:,ii) + eps) + eps);
    end
% Version = 2.7    
else
    % GantryRtn = KvSourceRtn - 90;
    angles = angles - 90;

    % Flip for interpolation
    if(Sec(2) - Sec(1) <0)
        Sec = flip(Sec);
        BlkAirNorm = flip(BlkAirNorm);
        gBlk = flip(gBlk, 3);
    end    
    
    % interpolation weights
    for ii = 1:length(angles)
        [left, weights] = interp_weight(angles(ii), Sec);
        interp_blk = weights(1) * gBlk(:,:,left) + weights(2) * gBlk(:,:,left+1);
        % Correction factor
        CF = airnorm(ii)/(0.5*BlkAirNorm(left) + 0.5*BlkAirNorm(left+1));
        gproj(:,:,ii) = log(CF*interp_blk./(gproj(:,:,ii)+eps) + eps);
    end
end

proj = gather(gproj);
% GPU clear
reset(gpuDevice(gpuids.devices(0)+1));

end
