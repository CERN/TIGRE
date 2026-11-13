function [image,residual,quality]=CGLS_operator(data,op,niter,varargin)
% Explicit opt-in CGLS; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('CGLS',data,op,niter,varargin{:});
end
