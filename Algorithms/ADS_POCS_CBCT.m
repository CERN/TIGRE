function [ f ] = ADS_POCS_CBCT(proj,geo,angles)
%ADS_POCS_CBCT Summary of this function goes here
%   Detailed explanation goes here

%% Create weigthing matrices

% Projection weigth, W
W=Ax(ones(geo.nVoxel'),geo,angles);  %
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;
% Back-Projection weigth, V
[x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
    -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
A = permute(angles, [1 3 2]);
V = (geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2;
V=sum(V,3);
clear A x y dx dz;
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;
%% vervatim implemetation of the paper
beta=1;
beta_red=0.995;
ng=20;
alpha=0.2;
rmax=0.95;
alpha_red=0.95;

f=zeros(geo.nVoxel');
g0=proj;
stop_criteria=0;
while ~stop_criteria %POCS
    f0=f;
    % 1 SART update (in the paper is an ART update)
    for jj=1:length(alpha);
        if size(offOrigin,2)==length(alpha)
            geo.OffOrigin=offOrigin(:,jj);
        end
         if size(offDetector,2)==length(alpha)
            geo.offDetector=offDetector(:,jj);
        end
        
        proj_err=proj(:,:,jj)-Ax(f,geo,alpha(jj));        %                                 (b-Ax)
        weighted_err=W(:,:,jj).*proj_err;                 %                          W^-1 * (b-Ax)
        backprj=Atb(weighted_err,geo,alpha(jj));          %                     At * W^-1 * (b-Ax)
        weigth_backprj=bsxfun(@times,1./V,backprj);       %                 V * At * W^-1 * (b-Ax)
        f=f+lambda*weigth_backprj;                        % x= x + lambda * V * At * W^-1 * (b-Ax) 
    end
    % Enforce positivity
    f(f<0)=0; 
    fres=f;
    g=Ax(f,geo,alpha(jj));
    dd=abs(g-g0);% ?
    dp=abs(f-f0);
    
    if iter==1
       dtvg=alpha*dp; 
    end
    f0=f;
    for ii=1:ng
        
    end
end



end

