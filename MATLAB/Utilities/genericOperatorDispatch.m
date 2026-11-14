function [handled,image,residual,quality]=genericOperatorDispatch(name,data,niter,args)
% Opt-in only. Calls without Operator retain their original implementation.
handled=false; image=[]; residual=[]; quality=[];
assert(mod(numel(args),2)==0,'TIGRE:operator:options','Use name/value options');
which=find(strcmpi(args(1:2:end),'Operator')); if isempty(which), return; end
assert(numel(which)==1,'TIGRE:operator:options','Supply one Operator');
k=2*which-1; op=args{k+1}; args(k:k+1)=[]; handled=true;
[image,residual,quality]=allOperatorRecon(name,data,op,niter,args{:});
end
