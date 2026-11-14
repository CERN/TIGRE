function [image,residual,quality]=MLEM_surface(data,geo,surface,niter,varargin)
% Explicit opt-in MLEM; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('MLEM',data,geo,surface,niter,varargin{:});
end
