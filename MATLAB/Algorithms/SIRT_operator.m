function [image,residual,quality]=SIRT_operator(data,op,niter,varargin)
% Explicit opt-in SIRT; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('SIRT',data,op,niter,varargin{:});
end
