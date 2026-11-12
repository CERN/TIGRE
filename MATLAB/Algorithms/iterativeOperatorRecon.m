function [x,residual]=iterativeOperatorRecon(data,op,niter,method,varargin)
% Shared CGLS/SIRT/OS-SART kernels; operators retain MATLAB array layouts.
p=inputParser; p.PartialMatching=false;
addParameter(p,'Init',[]); addParameter(p,'Relaxation',1); addParameter(p,'Nonnegative',false);
addParameter(p,'BlockSize',20); addParameter(p,'Tolerance',0);
addParameter(p,'TViter',20); addParameter(p,'TVLambda',50); addParameter(p,'RelaxationReduction',1);
parse(p,varargin{:}); a=p.Results;
validateattributes(niter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(a.Relaxation,{'numeric'},{'scalar','positive','finite'});
validateattributes(a.BlockSize,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(a.Tolerance,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(a.TViter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(a.TVLambda,{'numeric'},{'scalar','positive','finite'});
validateattributes(a.RelaxationReduction,{'numeric'},{'scalar','positive','finite'});
validateattributes(a.Nonnegative,{'logical'},{'scalar'});
b=operatorPrepareData(data,op);
x=zeros(op.volume_shape,'single');
if ~isempty(a.Init)
    assert(isnumeric(a.Init)&&isreal(a.Init)&&isequal([size(a.Init,1) size(a.Init,2) size(a.Init,3)],op.volume_shape),'TIGRE:operator:init','Invalid Init');
    x=single(a.Init); assert(all(isfinite(x(:))),'TIGRE:operator:init','Invalid Init');
end
residual=zeros(1,0); if niter==0, return; end
r=b-op.forward(x); initial=sqrt(squared(r));
if strcmp(method,'cgls')
    direction=op.adjoint(r); gamma=squared(direction);
    for k=1:niter
        if gamma==0 || (a.Tolerance>0 && sqrt(squared(r))<=a.Tolerance*initial), break; end
        q=op.forward(direction); qq=squared(q); if qq==0, break; end
        alpha=gamma/qq; x=x+single(alpha)*direction; r=r-single(alpha)*q;
        residual(end+1)=sqrt(squared(r)); next=op.adjoint(r); gnew=squared(next);
        direction=next+single(gnew/gamma)*direction; gamma=gnew;
    end
    return;
end
nv=op.projection_shape(3); block=a.BlockSize;
if strcmp(method,'sirt'), block=nv; end
nb=ceil(nv/block); subs=cell(1,nb); rows=subs; indices=subs;
for k=1:nb
    indices{k}=k:nb:nv; subs{k}=op.subset(indices{k});
    row=subs{k}.forward(ones(op.volume_shape,'single'));
    w=zeros(size(row),'single'); w(row>1e-6)=1./row(row>1e-6); rows{k}=w;
end
t=1; previous=x; relaxation=a.Relaxation;
for iteration=1:niter
    if a.Tolerance>0 && sqrt(squared(r))<=a.Tolerance*initial, break; end
    if strcmp(method,'fast_os_sart')
        tn=(1+sqrt(1+4*t*t))/2; saved=x; x=x+single((t-1)/tn)*(x-previous); previous=saved; t=tn;
    end
    for k=1:nb
        sub=subs{k}; valid=ones(sub.projection_shape,'single');
        if ~isempty(sub.valid_mask), valid=single(sub.valid_mask); end
        col=sub.adjoint(valid); v=zeros(size(col),'single'); v(col>1e-6)=1./col(col>1e-6);
        x=x+single(relaxation)*v.*sub.adjoint((b(:,:,indices{k})-sub.forward(x)).*rows{k});
        if a.Nonnegative, x=max(x,0); end
    end
    if contains(method,'_tv') && a.TViter>0 && max(x(:))>min(x(:))
        ids=GpuIds(); if isfield(op,'device'), ids.devices=int32(op.device); end
        x=im3DDenoise(x,'TV',a.TViter,a.TVLambda,'gpuids',ids);
        if a.Nonnegative, x=max(x,0); end
    end
    relaxation=relaxation*a.RelaxationReduction;
    r=b-op.forward(x); residual(end+1)=sqrt(squared(r));
end
end
function value=squared(x)
value=sum(double(x(:)).^2);
end
