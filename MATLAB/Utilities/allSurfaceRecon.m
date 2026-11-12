function [image,residual,quality]=allSurfaceRecon(name,data,geo,surface,niter,varargin)
args=varargin; found=find(strcmpi(args(1:2:end),'Device')); device=0;
assert(numel(found)<=1,'TIGRE:operator:options','Supply one Device');
if ~isempty(found), k=2*found-1; device=args{k+1}; args(k:k+1)=[]; end
op=projectionOperator(geo,surface,'Device',device);
[image,residual,quality]=allOperatorRecon(name,data,op,niter,args{:});
end
