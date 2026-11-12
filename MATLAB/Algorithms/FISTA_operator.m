function [image,residual,quality]=FISTA_operator(data,op,niter,varargin)
% Explicit opt-in FISTA; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('FISTA',data,op,niter,varargin{:});
end
