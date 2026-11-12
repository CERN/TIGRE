function [image,residual,quality]=OS_ASD_POCS_operator(data,op,niter,varargin)
% Explicit opt-in OS_ASD_POCS; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('OS_ASD_POCS',data,op,niter,varargin{:});
end
