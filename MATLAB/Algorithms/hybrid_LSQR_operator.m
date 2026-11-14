function [image,residual,quality]=hybrid_LSQR_operator(data,op,niter,varargin)
% Explicit opt-in hybrid_LSQR; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('hybrid_LSQR',data,op,niter,varargin{:});
end
