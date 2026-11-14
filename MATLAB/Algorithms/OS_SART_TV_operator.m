function [image,residual,quality]=OS_SART_TV_operator(data,op,niter,varargin)
% Explicit opt-in OS_SART_TV; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('OS_SART_TV',data,op,niter,varargin{:});
end
