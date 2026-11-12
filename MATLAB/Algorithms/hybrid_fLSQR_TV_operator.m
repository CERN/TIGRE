function [image,residual,quality]=hybrid_fLSQR_TV_operator(data,op,niter,varargin)
% Explicit opt-in hybrid_fLSQR_TV; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('hybrid_fLSQR_TV',data,op,niter,varargin{:});
end
