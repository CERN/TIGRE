function [image,residual,quality]=SART_TV_surface(data,geo,surface,niter,varargin)
% Explicit opt-in SART_TV; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('SART_TV',data,geo,surface,niter,varargin{:});
end
