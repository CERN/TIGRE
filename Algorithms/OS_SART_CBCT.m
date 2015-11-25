function [res,errorL2]=OS_SART_CBCT(proj,geo,alpha,niter,block_size,lambda)

%% Deal with input parameters
if nargin<6
    lambda=1;
end

%% Create weigthing matrices

% Projection weigth, W
% Projection weigth, W
W=Ax(ones(geo.nVoxel'),geo,alpha);  %
W(W<min(geo.dVoxel))=Inf;
W=1./W;
% Back-Projection weigth, V
[x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
    -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
A = permute(alpha, [1 3 2]);
V = (geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2;
V=sum(V,3);
clear A x y dx dz;


%% initialize image With FDK

% TODO : this should be optional
% TODO : Actually call FDK here
% geo.filter='ram-lak';
% proj_filt = filtering(proj,geo,alpha);
% geo=rmfield(geo,'filter');
% res=Atb(proj_filt,geo,alpha);
res=zeros(geo.nVoxel');


%% Iterate
errorL2=norm(proj(:));
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;

% TODO : Add options for Stopping criteria
for ii=1:niter
    if ii==1;tic;end
    
    for jj=1:block_size:length(alpha);
        % idex of the Oriented subsets
        range=jj:block_size+jj-1;
        range(range>length(alpha))=[]; % for the last subset
        
        if size(offOrigin,2)==length(alpha)
            geo.OffOrigin=offOrigin(:,range);
        end
         if size(offDetector,2)==length(alpha)
            geo.offDetector=offDetector(:,range);
        end
        
        proj_err=proj(:,:,range)-Ax(res,geo,alpha(range));      %                                 (b-Ax)
        weighted_err=W(:,:,range).*proj_err;                    %                          W^-1 * (b-Ax)
        backprj=Atb(weighted_err,geo,alpha(range));             %                     At * W^-1 * (b-Ax)
        weigth_backprj=bsxfun(@times,1./V,backprj);             %                 V * At * W^-1 * (b-Ax)
        res=res+lambda*weigth_backprj;                          % x= x + lambda * V * At * W^-1 * (b-Ax)
        
        
        
    end
    errornow=norm(proj_err(:));                           % Compute error norm2 of b-Ax
    % If the error is not minimized (Armijo rules)
    if errornow>errorL2(end)
        return;
    end
    errorL2=[errorL2 errornow];
    
    if ii==1;
        expected_time=toc*(niter-1);   
        disp('OS-SART');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
end





end