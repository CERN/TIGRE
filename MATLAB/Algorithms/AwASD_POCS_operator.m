function [image,residual,quality]=AwASD_POCS_operator(data,op,niter,varargin)
% Explicit opt-in AwASD_POCS; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('AwASD_POCS',data,op,niter,varargin{:});
end
