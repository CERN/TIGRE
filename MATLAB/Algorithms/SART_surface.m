function [image,residual,quality]=SART_surface(data,geo,surface,niter,varargin)
% Explicit opt-in SART; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('SART',data,geo,surface,niter,varargin{:});
end
