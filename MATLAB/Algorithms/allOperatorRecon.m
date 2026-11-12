function [image,residual,quality]=allOperatorRecon(name,data,op,niter,varargin)
% Shared public operator dispatcher. Conventional methods opt in explicitly.
assert(ismember(lower(name),lower({'CGLS', 'SIRT', 'SART', 'OS_SART', 'fast_OS_SART', 'SART_TV', 'OS_SART_TV', 'iterativereconalg', 'LSQR', 'LSMR', 'hybrid_LSQR', 'IRN_TV_CGLS', 'hybrid_fLSQR_TV', 'AB_GMRES', 'BA_GMRES', 'FISTA', 'ISTA', 'MLEM', 'ASD_POCS', 'OS_ASD_POCS', 'AwASD_POCS', 'OS_AwASD_POCS', 'PCSD', 'OS_PCSD', 'AwPCSD', 'OS_AwPCSD'})), 'TIGRE:operator:unsupported','Method is not included in this operator slice');
quality=[]; key=lower(name); args=varargin;
pocs={'asd_pocs','os_asd_pocs','awasd_pocs','os_awasd_pocs','pcsd','os_pcsd','awpcsd','os_awpcsd'};
if ismember(key,pocs)
    weighted=contains(key,'aw'); ordered=startsWith(key,'os_'); pcsd=contains(key,'pcsd');
    [~,image,residual]=operatorPocsDispatch(data,niter,weighted,ordered,{'Operator',op,args{:}},pcsd);
elseif ismember(key,{'sart','sart_tv','os_sart_tv','iterativereconalg'})
    for k=1:2:numel(args)
        switch lower(args{k})
            case 'lambda', args{k}='Relaxation';
            case 'lambda_red', args{k}='RelaxationReduction';
            case 'nonneg', args{k}='Nonnegative';
        end
    end
    if ~any(strcmpi(args(1:2:end),'Nonnegative')), args=[args {'Nonnegative',true}]; end
    if ismember(key,{'sart','sart_tv'}), args=[args {'BlockSize',1}]; end
    [image,residual]=iterativeOperatorRecon(data,op,niter,key,args{:});
elseif ismember(key,{'cgls','sirt','os_sart','fast_os_sart'})
    [image,residual]=iterativeOperatorRecon(data,op,niter,key,args{:});
elseif ismember(key,{'fista','ista','mlem'})
    [image,residual]=proximalOperatorRecon(data,op,niter,key,args{:});
else
    [image,residual,quality]=existingOperatorRecon(name,data,op,niter,args{:});
end
end
