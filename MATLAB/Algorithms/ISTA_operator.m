function [image,residual,quality]=ISTA_operator(data,op,niter,varargin)
% Explicit opt-in ISTA; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('ISTA',data,op,niter,varargin{:});
end
