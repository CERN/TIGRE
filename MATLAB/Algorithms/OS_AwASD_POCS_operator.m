function [image,residual,quality]=OS_AwASD_POCS_operator(data,op,niter,varargin)
% Explicit opt-in OS_AwASD_POCS; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('OS_AwASD_POCS',data,op,niter,varargin{:});
end
