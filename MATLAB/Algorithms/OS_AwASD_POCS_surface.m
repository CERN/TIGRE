function [image,residual,quality]=OS_AwASD_POCS_surface(data,geo,surface,niter,varargin)
% Explicit opt-in OS_AwASD_POCS; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('OS_AwASD_POCS',data,geo,surface,niter,varargin{:});
end
