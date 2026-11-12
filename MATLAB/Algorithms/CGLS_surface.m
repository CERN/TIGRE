function [image,residual,quality]=CGLS_surface(data,geo,surface,niter,varargin)
% Explicit opt-in CGLS; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('CGLS',data,geo,surface,niter,varargin{:});
end
