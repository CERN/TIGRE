function [image,residual,quality]=SIRT_surface(data,geo,surface,niter,varargin)
% Explicit opt-in SIRT; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('SIRT',data,geo,surface,niter,varargin{:});
end
