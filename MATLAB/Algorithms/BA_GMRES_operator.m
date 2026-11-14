function [image,residual,quality]=BA_GMRES_operator(data,op,niter,varargin)
% Explicit opt-in BA_GMRES; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('BA_GMRES',data,op,niter,varargin{:});
end
