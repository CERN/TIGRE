function [image,residual,quality]=OS_AwPCSD_operator(data,op,niter,varargin)
% Explicit opt-in OS_AwPCSD; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('OS_AwPCSD',data,op,niter,varargin{:});
end
