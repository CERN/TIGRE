function [image,residual,quality]=IRN_TV_CGLS_operator(data,op,niter,varargin)
% Explicit opt-in IRN_TV_CGLS; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('IRN_TV_CGLS',data,op,niter,varargin{:});
end
