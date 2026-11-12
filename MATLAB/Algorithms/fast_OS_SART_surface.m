function [image,residual,quality]=fast_OS_SART_surface(data,geo,surface,niter,varargin)
% Explicit opt-in fast_OS_SART; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('fast_OS_SART',data,geo,surface,niter,varargin{:});
end
