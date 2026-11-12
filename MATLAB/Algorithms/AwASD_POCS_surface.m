function [image,residual,quality]=AwASD_POCS_surface(data,geo,surface,niter,varargin)
% Explicit opt-in AwASD_POCS; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('AwASD_POCS',data,geo,surface,niter,varargin{:});
end
