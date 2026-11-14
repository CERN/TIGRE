function [image,residual,quality]=MLEM_operator(data,op,niter,varargin)
% Explicit opt-in MLEM; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('MLEM',data,op,niter,varargin{:});
end
