function [image,residual,quality]=OS_AwPCSD_surface(data,geo,surface,niter,varargin)
% Explicit opt-in OS_AwPCSD; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('OS_AwPCSD',data,geo,surface,niter,varargin{:});
end
