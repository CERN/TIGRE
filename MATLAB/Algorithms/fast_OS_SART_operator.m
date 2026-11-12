function [image,residual,quality]=fast_OS_SART_operator(data,op,niter,varargin)
% Explicit opt-in fast_OS_SART; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('fast_OS_SART',data,op,niter,varargin{:});
end
