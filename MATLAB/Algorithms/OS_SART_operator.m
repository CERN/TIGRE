function [image,residual,quality]=OS_SART_operator(data,op,niter,varargin)
% Explicit opt-in OS_SART; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('OS_SART',data,op,niter,varargin{:});
end
