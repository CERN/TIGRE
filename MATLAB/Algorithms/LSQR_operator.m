function [image,residual,quality]=LSQR_operator(data,op,niter,varargin)
% Explicit opt-in LSQR; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('LSQR',data,op,niter,varargin{:});
end
