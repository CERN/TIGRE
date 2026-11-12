function [image,residual,quality]=OS_SART_TV_surface(data,geo,surface,niter,varargin)
% Explicit opt-in OS_SART_TV; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('OS_SART_TV',data,geo,surface,niter,varargin{:});
end
