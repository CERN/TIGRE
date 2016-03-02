function [res,errorL2,qualMeas]=OS_SART_CBCT(proj,geo,alpha,niter,varargin)
% OS_SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets
%              Simultaneous Algebraic Reconxtruction Techique algorithm
%
%   OS_SART_CBCT(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
%
%   OS_SART_CBCT(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%   'BlockSize':   Sets the projection block size used simultaneously. If
%                  BlockSize = 1 OS-SART becomes SART and if  BlockSize = length(alpha)
%                  then OS-SART becomes SIRT. Default is alpha/20.
%
%   'lambda':      Sets the value of the hyperparameter. Default is 1
%
%   'lambdared':   Reduction of lambda.Every iteration
%                  lambda=lambdared*lambda. Default is 0.95
%
%   'Init':        Describes diferent initialization techniques.
%                  'none'     : Initializes the image to zeros (default)
%                  'FDK'      : intializes image to FDK reconstrucition
%                  'multigrid': Initializes image by solving the problem in
%                               small scale and increasing it when relative
%                               convergence is reached.
%                  'image'    : Initialization using a user specified
%                               image. Not recomended unless you really
%                               know what you are doing.
%   'InitImg'      an image for the 'image' initialization. Aviod.
%
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
%   'QualMeas'     Asks the algorithm for a set of quality measurement
%                  parameters. Input should contain a cell array of desired
%                  quality measurement names. Example: {'CC','RMSE','MSSIM'}
%                  These will be computed in each iteration. 
%
%   
% OUTPUTS:
%    [img]                       will output the reconstructed image
%    [img,errorL2]               will output the L2 norm of the residual 
%                                (the function being minimized)
%    [img,errorL2,qualMease]     will output the quality measurements asked
%                                by the input 'QualMeas'
%
%
%% Deal with input parameters

opts=     {'BlockSize','lambda','Init','InitImg','Verbose','lambdaRed'};
defaults= [   1  ,  1  ,    1   ,1 ,1,1];



% Check inputs
nVarargs = length(varargin);
if mod(nVarargs,2)
    error('CBCT:plotImgs:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,varargin{ii}));
    if ~isempty(ind)
        defaults(ind)=0;
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extranc value from input
   if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,varargin{jj}));
            jj=jj+1;
        end
        val=varargin{jj};
    end
    
    switch opt
        % % % % % % % Verbose
        case 'Verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
        % % % % % % % hyperparameter, LAMBDA
        case 'lambda'
            if default
                lambda=0.95;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid lambda')
                end
                lambda=val;
            end
         case 'lambdaRed'
            if default
                lamdbared=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid lambda')
                end
                lamdbared=val;
            end
        case 'BlockSize'
            if default
                block_size=20;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid BlockSize')
                end
                block_size=val;
            end
            
        case 'Init'
            res=[];
            if default || strcmp(val,'none')
                res=zeros(geo.nVoxel');
                continue;
            end
            if strcmp(val,'FDK')
                res=FDK_CBCT(proj,geo,alpha);
                continue;
            end
            if strcmp(val,'multigrid')
                res=init_multigrid(proj,geo,alpha);
                continue;
            end
            if strcmp(val,'image')
                initwithimage=1;
                continue;
            end
            if isempty(res)
               error('CBCT:OS_SART_CBCT:InvalidInput','Invalid Init option') 
            end
            % % % % % % % ERROR
        case 'InitImg'
            if default
                continue;
            end
            if exist('initwithimage','var');
                if isequal(size(val),geo.nVoxel');
                    res=val;
                else
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid image for initialization');
                end
            end
        otherwise
            error('CBCT:OS_SART_CBCT:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in OS_SART_CBCT()']);
    end
end









%% Create weighting matrices


% Projection weight, W
W=Ax(ones(geo.nVoxel'),geo,alpha,'Krylov');  
% What is Krylov?
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;

% Back-Projection weight, V
[x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
    -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
A = permute(alpha, [1 3 2]);
V = (geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2;
%C = bsxfun(fun,A,B) applies the element-by-element binary operation specified by the function handle fun to arrays A and B, with singleton expansion enabled. 
%@times is Array multiply
V=sum(V,3);
clear A x y dx dz;




%% Iterate
errorL2=norm(proj(:)); %Calculate for 2-norm
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;



% TODO : Add options for Stopping criteria
for ii=1:niter
    if (ii==1 && verbose==1);tic;end
    
    for jj=1:block_size:length(alpha);
        % index of the Oriented subsets

        range=jj:block_size+jj-1;
        
        range(range>length(alpha))=[]; % for the last subset
        
        if size(offOrigin,2)==length(alpha)
            geo.offOrigin=offOrigin(:,range);
        end
        if size(offDetector,2)==length(alpha)
            geo.offDetector=offDetector(:,range);
        end
        
        %proj is data: b=Ax
        %res= initial image is zero (default)
        proj_err=proj(:,:,range)-Ax(res,geo,alpha(range),'Krylov');      %                                 (b-Ax)
        weighted_err=W(:,:,range).*proj_err;                             %                          W^-1 * (b-Ax)
        backprj=Atb(weighted_err,geo,alpha(range));                      %                     At * W^-1 * (b-Ax)
        weigth_backprj=bsxfun(@times,1./V,backprj);                      %                 V * At * W^-1 * (b-Ax)
        
        
        rm=RMSE(res,res+lambda*weigth_backprj);   
        corr=CC(res,res+lambda*weigth_backprj);
        mss=MSSIM(res,res+lambda*weigth_backprj);
        res=res+lambda*weigth_backprj;                                   % x= x + lambda * V * At * W^-1 * (b-Ax)
        
        % Non-negativity constrain
        res(res<0)=0;

        %Store the values every iteration
        rmtotal(ii)=[rm];
        corrtotal(ii)=[corr];
        msstotal(ii)=[mss];
        
        
    end
    lambda=lambda*lamdbared;
    errornow=norm(proj_err(:));                           % Compute error norm2 of b-Ax
%     If the error is not minimized 
    if errornow>errorL2(end)*1.1
        return;
    end
    errorL2=[errorL2 errornow];
    
    if ii==1 && verbose==1
        
        expected_time=toc*(niter-1);
        disp('OS-SART');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    
 
    
 
end

 
 



end
function initres=init_multigrid(proj,geo,alpha)

finalsize=geo.nVoxel;
% start with 64
geo.nVoxel=[64;64;64];
geo.dVoxel=geo.sVoxel./geo.nVoxel;
if any(finalsize<geo.nVoxel)
    initres=zeros(finalsize');
    return;
end
niter=100;
nblock=20;
initres=zeros(geo.nVoxel');
while ~isequal(geo.nVoxel,finalsize)
    
    
    % solve subsampled grid
    initres=OS_SART_CBCT(proj,geo,alpha,niter,'BlockSize',nblock,'Init','image','InitImg',initres,'Verbose',0);
    
    % Get new dims.
    geo.nVoxel=geo.nVoxel*2;
    geo.nVoxel(geo.nVoxel>finalsize)=finalsize(geo.nVoxel>finalsize);
    geo.dVoxel=geo.sVoxel./geo.nVoxel;
    % Upsample!
    % (hopefully computer has enough memory............)
    [y, x, z]=ndgrid(linspace(1,size(initres,1),geo.nVoxel(1)),...
                     linspace(1,size(initres,2),geo.nVoxel(2)),...
                     linspace(1,size(initres,3),geo.nVoxel(3)));
    initres=interp3(initres,x,y,z);
    clear x y z 
end


end
