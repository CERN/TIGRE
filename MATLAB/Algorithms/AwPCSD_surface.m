function [image,residual,quality]=AwPCSD_surface(data,geo,surface,niter,varargin)
% Explicit opt-in AwPCSD; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('AwPCSD',data,geo,surface,niter,varargin{:});
end
