function [image,residual,quality]=SART_operator(data,op,niter,varargin)
% Explicit opt-in SART; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('SART',data,op,niter,varargin{:});
end
