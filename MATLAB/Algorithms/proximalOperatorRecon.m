function [x,history]=proximalOperatorRecon(data,op,niter,method,varargin)
% Operator-aware MLEM and ISTA/FISTA volume-TV steps.
p=inputParser; p.PartialMatching=false;
addParameter(p,'Init',[]); addParameter(p,'Hyper',[]); addParameter(p,'TViter',20);
addParameter(p,'lambda',.1); addParameter(p,'Nonneg',true);
addParameter(p,'fista_p',1); addParameter(p,'fista_q',1);
parse(p,varargin{:}); a=p.Results;
validateattributes(niter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(a.TViter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(a.lambda,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(a.Nonneg,{'logical'},{'scalar'});
validateattributes(a.fista_p,{'numeric'},{'scalar','positive','finite'});
validateattributes(a.fista_q,{'numeric'},{'scalar','nonnegative','finite'});
if ~isempty(a.Hyper), validateattributes(a.Hyper,{'numeric'},{'scalar','positive','finite'}); end
b=operatorPrepareData(data,op);
if strcmp(method,'mlem')
    assert(all(b(:)>=0),'TIGRE:operator:mlem','MLEM needs nonnegative data');
end
x=zeros(op.volume_shape,'single'); if strcmp(method,'mlem'), x=x+1; end
if ~isempty(a.Init)
    assert(isnumeric(a.Init)&&isreal(a.Init)&&isequal([size(a.Init,1) size(a.Init,2) size(a.Init,3)],op.volume_shape),'TIGRE:operator:init','Invalid Init');
    x=single(a.Init); assert(all(isfinite(x(:))),'TIGRE:operator:init','Invalid Init');
end
if a.Nonneg || strcmp(method,'mlem'), x=max(x,0); end
history=zeros(1,0); if niter==0, return; end
if strcmp(method,'mlem')
    assert(all(b(:)>=0),'TIGRE:operator:mlem','MLEM needs nonnegative data');
    valid=ones(op.projection_shape,'single'); if ~isempty(op.valid_mask), valid=single(op.valid_mask); end
    sensitivity=op.adjoint(valid); sensitivity(sensitivity<=0)=inf;
    for k=1:niter
        den=op.forward(x); den(den<=1e-6)=inf;
        x=max(x.*(op.adjoint(min(max(b./den,0),10))./sensitivity),0);
        r=op.forward(x)-b; history(end+1)=norm(double(r(:)));
    end
    return;
end
L=a.Hyper;
if isempty(L)
    saved=rng; cleanup=onCleanup(@() rng(saved)); rng(0);
    v=single(randn(op.volume_shape)); v=v/single(norm(double(v(:))));
    for k=1:8
        v=op.adjoint(op.forward(v)); scale=norm(double(v(:)));
        if scale==0, break; end
        v=v/single(scale);
    end
    L=max(2.1*scale,eps);
end
validateattributes(L,{'numeric'},{'scalar','positive','finite'});
ids=GpuIds(); if isfield(op,'device'), ids.devices=int32(op.device); end
prior=x; t=1;
for k=1:niter
    x=x+single(2/L)*op.adjoint(b-op.forward(x)); if a.Nonneg, x=max(x,0); end
    denoised=x;
    if a.lambda>0 && a.TViter>0 && max(x(:))>min(x(:))
        denoised=im3DDenoise(x,'TV',a.TViter,L/(2*a.lambda),'gpuids',ids);
    end
    if strcmp(method,'fista')
        next=(a.fista_p+sqrt(a.fista_q+4*t*t))/2;
        x=denoised+single((t-1)/next)*(denoised-prior); t=next; prior=denoised;
    else, x=denoised; end
    if a.Nonneg, x=max(x,0); end
    r=op.forward(x)-b; history(end+1)=norm(double(r(:)));
end
end
