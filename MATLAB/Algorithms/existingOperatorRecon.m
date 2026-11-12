function [image,residual,quality]=existingOperatorRecon(name,data,op,niter,varargin)
% Reuse original Krylov/statistical recurrences with bound operator projectors.
validateattributes(niter,{'numeric'},{'scalar','integer','nonnegative','finite'});
b=operatorPrepareData(data,op);
args=varargin;
assert(mod(numel(args),2)==0,'TIGRE:operator:options','Use name/value options');
keys=lower(string(args(1:2:end))); initial=[];
allowed={'init','verbose','qualmeas','gpuids','groundtruth'};
if ismember(lower(name),{'lsqr','lsmr','irn_tv_cgls'}), allowed=[allowed {'restart'}]; end
if ismember(lower(name),{'hybrid_lsqr','hybrid_flsqr_tv','irn_tv_cgls','lsmr'}), allowed=[allowed {'lambda'}]; end
if ismember(lower(name),{'hybrid_lsqr','hybrid_flsqr_tv'}), allowed=[allowed {'noiselevel'}]; end
if strcmpi(name,'IRN_TV_CGLS'), allowed=[allowed {'niter_outer'}]; end
if ismember(lower(name),{'ab_gmres','ba_gmres'}), allowed=[allowed {'backprojector'}]; end
assert(all(ismember(keys,string(allowed))),'TIGRE:operator:options','Unsupported operator solver option');
if ismember(lower(name),{'hybrid_lsqr','hybrid_flsqr_tv','irn_tv_cgls','lsmr'}) && ~any(keys=='lambda')
    args=[args {'lambda',1}];
end
if strcmpi(name,'IRN_TV_CGLS')
    entry=find(keys=='niter_outer'); outer=4;
    if isempty(entry), args=[args {'niter_outer',outer}]; else, outer=args{2*entry}; end
    validateattributes(outer,{'numeric'},{'scalar','integer','positive','finite'});
    niter=niter*outer;
end
for k=1:numel(keys)
    value=args{2*k};
    if keys(k)=='init' && (ischar(value)||isstring(value)) && ~strcmpi(value,'none')
        error('TIGRE:operator:init','FDK/multigrid initializers are geometry-specific');
    end
    if keys(k)=='backprojector' && ~strcmpi(value,'matched')
        error('TIGRE:operator:backprojector','Use the operator adjoint, not FDK');
    end
    if keys(k)=='init' && isnumeric(value)
        assert(isequal([size(value,1) size(value,2) size(value,3)],op.volume_shape),'TIGRE:operator:init','Invalid Init');
        initial=single(value); args{2*k}='image'; args=[args {'InitImg',initial}];
        assert(isreal(initial)&&all(isfinite(initial(:))),'TIGRE:operator:init','Invalid finite Init');
    end
    if ismember(keys(k),{'lambda','noiselevel'})
        validateattributes(value,{'numeric'},{'scalar','real','nonnegative','finite'});
    end
end
if ~any(keys=='verbose'), args=[args {'Verbose',0}]; end
if ~any(keys=='gpuids')
    ids=GpuIds(); if isfield(op,'device'), ids.devices=int32(op.device); end
    args=[args {'gpuids',ids}];
end
geo=struct('nVoxel',double(op.volume_shape(:)),'nDetector',double(op.projection_shape(1:2))', ...
           'operator_binding',op);
image=zeros(op.volume_shape,'single'); residual=zeros(1,0); quality=[];
if ~isempty(initial), image=initial; end
if niter==0 || (~any(b(:)) && ~any(keys=='init')), return; end
if ismember(lower(name),{'lsqr','lsmr','hybrid_lsqr','ab_gmres','ba_gmres','hybrid_flsqr_tv'})
    r=b; if ~isempty(initial), r=b-op.forward(initial); end
    stationary=~any(r(:));
    if ~stationary, normal=op.adjoint(r); stationary=~any(normal(:)); end
    if stationary, return; end
end
angles=1:op.projection_shape(3);
if nargout(name)>=3
    [image,residual,quality]=feval(name,b,geo,angles,niter,args{:});
else
    [image,quality]=feval(name,b,geo,angles,niter,args{:});
    r=op.forward(single(image))-b; residual=norm(double(r(:)));
end
image=single(image);
assert(all(isfinite(image(:))),'TIGRE:operator:finite','Algorithm produced nonfinite values');
end
