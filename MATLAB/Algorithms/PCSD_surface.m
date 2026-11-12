function [image,residual,quality]=PCSD_surface(data,geo,surface,niter,varargin)
% Explicit opt-in PCSD; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('PCSD',data,geo,surface,niter,varargin{:});
end
