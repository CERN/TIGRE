function [image,residual,quality]=OS_ASD_POCS_surface(data,geo,surface,niter,varargin)
% Explicit opt-in OS_ASD_POCS; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('OS_ASD_POCS',data,geo,surface,niter,varargin{:});
end
