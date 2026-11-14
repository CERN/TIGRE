function [image,residual,quality]=LSMR_surface(data,geo,surface,niter,varargin)
% Explicit opt-in LSMR; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('LSMR',data,geo,surface,niter,varargin{:});
end
