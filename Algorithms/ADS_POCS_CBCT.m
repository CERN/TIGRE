function [ fres,tv ] = ADS_POCS_CBCT(proj,geo,angles,maxiter,epsilon)
%ADS_POCS_CBCT Summary of this function goes here
%   Detailed explanation goes here
verbose=1;
%% Create weigthing matrices

% Projection weigth, W
W=Ax(ones(geo.nVoxel'),geo,angles,'Krylov');  %
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
block_size=20;
%% vervatim implemetation of the paper
beta=1;
beta_red=0.995;
ng=50;
alpha=0.005;
rmax=0.95;
alpha_red=0.95;

f=zeros(geo.nVoxel');
g0=proj;
stop_criteria=0;
iter=0;
while ~stop_criteria %POCS
    if ~iter;tic; end
    iter=iter+1;
    f0=f;
    % 1 OS-SART update (in the paper is an ART update)
        for jj=1:block_size:length(angles);
            % idex of the Oriented subsets
            range=jj:block_size+jj-1;
            range(range>length(angles))=[]; % for the last subset

            if size(offOrigin,2)==length(angles)
                geo.offOrigin=offOrigin(:,range);
            end
            if size(offDetector,2)==length(angles)
                geo.offDetector=offDetector(:,range);
            end

            proj_err=proj(:,:,range)-Ax(f,geo,angles(range),'Krylov');      %                                 (b-Ax)
            weighted_err=W(:,:,range).*proj_err;                            %                          W^-1 * (b-Ax)
            backprj=Atb(weighted_err,geo,angles(range));                    %                     At * W^-1 * (b-Ax)
            weigth_backprj=bsxfun(@times,1./V,backprj);                     %                 V * At * W^-1 * (b-Ax)
            f=f+beta*weigth_backprj;                                        % x= x + lambda * V * At * W^-1 * (b-Ax)

        end
    geo.offDetector=offDetector;
    geo.offOrigin=offOrigin;
    % Enforce positivity
    f(f<0)=0; 
    fres=f;
    g=Ax(f,geo,angles);
    dd=im3Dnorm(g-g0,'L2');
    dp_vec=(f-f0);
    dp=im3Dnorm(dp_vec,'L2');
    
    if iter==1
       dtvg=alpha*dp; 
    end
    f0=f;
%%     This is the MATLAB CODE, the functions are sill in the library, but CUDA is used nowadays
%     for ii=1:ng
%         % Steepest descend of TV norm
%         tv(ng*(iter-1)+ii)=im3Dnorm(f,'TV','forward');
%         df=gradientTVnorm(f,'forward');
%         df=df./im3Dnorm(df,'L2');
%         f=f-dtvg.*df;    
%     end
%%
    f=minimizeTV(f0,dtvg,ng);
    dg_vec=(f-f0);
    dg=im3Dnorm(dg_vec,'L2');
    if dg>rmax*dp &&dd>epsilon
        dtvg=dtvg*alpha_red;
    end
   beta=beta*beta_red;
   %check stop criteria
   c=dot(dg_vec(:),dp_vec(:))/(norm(dg_vec(:),2)*norm(dp_vec(:),2));
   if (c<-0.99 && dd<=epsilon) || beta<0.005|| iter>maxiter
       stop_criteria=true;
   end
   if (ii==1 && verbose==1);
        expected_time=toc*niter;   
        disp('ADS-POCS');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end

end



end

