function [x,history]=ASD_POCS_operator(data,op,niter,varargin)
%ASD_POCS_OPERATOR Shared ASD/OS/weighted POCS implementation.
% Required MaxL2err: absolute masked L2 tolerance (no FDK-derived default).
% Options mirror Python: lambda, lambda_red, alpha, alpha_red, Ratio,
% TViter, Nonneg, Init. BlockSize=1, Weighted=false, Delta=-.005.
p=inputParser; p.PartialMatching=false;
addParameter(p,'MaxL2err',[]); addParameter(p,'lambda',1); addParameter(p,'lambda_red',.99);
addParameter(p,'alpha',.002); addParameter(p,'alpha_red',.95); addParameter(p,'Ratio',.95);
addParameter(p,'TViter',20); addParameter(p,'Nonneg',true); addParameter(p,'Init',[]);
addParameter(p,'BlockSize',1); addParameter(p,'Weighted',false); addParameter(p,'Delta',-.005);
addParameter(p,'PCSD',false);
parse(p,varargin{:}); a=p.Results;
validateattributes(niter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(a.MaxL2err,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(a.BlockSize,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(a.TViter,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(a.lambda,{'numeric'},{'scalar','positive','finite'});
for name={'lambda_red','alpha_red'}
    validateattributes(a.(name{1}),{'numeric'},{'scalar','>=',0,'<=',1,'finite'});
end
validateattributes(a.alpha,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(a.Ratio,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(a.Nonneg,{'logical'},{'scalar'}); validateattributes(a.Weighted,{'logical'},{'scalar'});
validateattributes(a.PCSD,{'logical'},{'scalar'});
if a.Weighted
    validateattributes(a.Delta,{'numeric'},{'real','vector','finite','nonzero'});
    assert(isscalar(a.Delta),'TIGRE:operator:delta','Directional Delta requires a separate native TV extension; use scalar Delta');
    delta=single(a.Delta);
    assert(all(isfinite(delta(:)))&&all(delta(:)~=0),'TIGRE:operator:delta','Invalid float32 Delta');
end
b=operatorPrepareData(data,op); shape=op.projection_shape;
x=zeros(op.volume_shape,'single');
if ~isempty(a.Init)
    assert(isnumeric(a.Init)&&isreal(a.Init)&&isequal([size(a.Init,1) size(a.Init,2) size(a.Init,3)],op.volume_shape),'TIGRE:operator:init','Init must be a volume');
    x=single(a.Init); assert(all(isfinite(x(:))),'TIGRE:operator:init','Invalid Init');
end
if a.Nonneg, x=max(x,0); end
keys={'residual','data_residual','beta','tv_step','cosine'};
for k=1:numel(keys), history.(keys{k})=zeros(1,0); end
if niter==0, return; end
nblocks=ceil(shape(3)/a.BlockSize); blocks=cell(1,nblocks); subs=blocks; weights=blocks;
for k=1:nblocks
    blocks{k}=k:nblocks:shape(3); subs{k}=op.subset(blocks{k});
    row=subs{k}.forward(ones(op.volume_shape,'single'));
    w=zeros(size(row),'single'); w(row>1e-6)=1./row(row>1e-6); weights{k}=w;
end
ids=[];
beta=a.lambda; step=0; referenceResidual=[];
for iteration=1:niter
    previous=x;
    active=1:nblocks;
    if a.PCSD
        before=length2(op.forward(x)-b);
        if before<=a.MaxL2err, active=[]; end
    end
    for k=active
        sub=subs{k}; valid=ones(sub.projection_shape,'single');
        if ~isempty(sub.valid_mask), valid=single(sub.valid_mask); end
        col=sub.adjoint(valid); v=zeros(size(col),'single'); v(col>1e-6)=1./col(col>1e-6);
        r=(b(:,:,blocks{k})-sub.forward(x)).*weights{k};
        x=x+single(beta)*v.*sub.adjoint(r); if a.Nonneg, x=max(x,0); end
    end
    dd=length2(op.forward(x)-b); dpvec=x-previous; dp=length2(dpvec);
    if a.PCSD
        if iteration==1
            step=double(before>0); referenceResidual=dd;
        elseif referenceResidual>0, step=before/referenceResidual;
        else, step=0; end
    elseif iteration==1, step=a.alpha*dp; end
    previous=x; used=step;
    if a.TViter>0 && step>0
        if isempty(ids)
            ids=GpuIds(); if isfield(op,'device'), ids.devices=int32(op.device); end
        end
        if a.Weighted, x=minimizeAwTV(x,step,a.TViter,delta,'gpuids',ids);
        else, x=minimizeTV(x,step,a.TViter,'gpuids',ids); end
        if a.Nonneg, x=max(x,0); end
    end
    dgvec=x-previous; dg=length2(dgvec);
    if ~a.PCSD && dg>a.Ratio*dp && dd>a.MaxL2err, step=step*a.alpha_red; end
    cosine=sum(double(dgvec(:)).*double(dpvec(:)))/max(dg*dp,1e-6);
    history.residual(end+1)=length2(op.forward(x)-b); history.data_residual(end+1)=dd;
    history.beta(end+1)=beta; history.tv_step(end+1)=used; history.cosine(end+1)=cosine;
    beta=beta*a.lambda_red;
    if (cosine<-.99 && dd<=a.MaxL2err)||beta<.005, break; end
end
end
function n=length2(x)
n=sqrt(sum(double(x(:)).^2));
end
