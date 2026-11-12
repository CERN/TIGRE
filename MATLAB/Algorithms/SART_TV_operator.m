function [image,residual,quality]=SART_TV_operator(data,op,niter,varargin)
% Explicit opt-in SART_TV; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('SART_TV',data,op,niter,varargin{:});
end
