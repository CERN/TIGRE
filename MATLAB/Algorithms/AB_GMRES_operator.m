function [image,residual,quality]=AB_GMRES_operator(data,op,niter,varargin)
% Explicit opt-in AB_GMRES; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('AB_GMRES',data,op,niter,varargin{:});
end
