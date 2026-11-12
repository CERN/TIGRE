function [handled,image,history]=operatorPocsDispatch(data,niter,weighted,ordered,args,pcsd)
% Internal opt-in bridge; no Operator option means the legacy path runs.
handled=false; image=[]; history=[];
if nargin<6, pcsd=false; end
assert(mod(numel(args),2)==0,'TIGRE:operator:options','Use name/value options');
which=find(strcmpi(args(1:2:end),'Operator'));
if isempty(which), return; end
assert(numel(which)==1,'TIGRE:operator:options','Supply one Operator');
k=2*which-1; op=args{k+1}; args(k:k+1)=[]; handled=true;
if pcsd
    allowed={'maxl2err','lambda','lambda_red','tviter','nonneg','init'};
    if ordered, allowed=[allowed {'blocksize'}]; end
    if weighted, allowed=[allowed {'delta'}]; end
    assert(all(ismember(lower(string(args(1:2:end))),string(allowed))), ...
        'TIGRE:operator:options','Unsupported PCSD operator option');
end
if ordered, block=20; else, block=1; end
if ~ordered
    entry=find(strcmpi(args(1:2:end),'BlockSize'));
    for item=entry
        assert(isequal(args{2*item},1),'TIGRE:operator:block','Use the OS variant for BlockSize > 1');
    end
end
[image,history]=ASD_POCS_operator(data,op,niter,'BlockSize',block,args{:},'Weighted',weighted,'PCSD',pcsd);
end
