function [image,residual,quality]=iterativereconalg_surface(data,geo,surface,niter,varargin)
% Explicit opt-in iterativereconalg; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('iterativereconalg',data,geo,surface,niter,varargin{:});
end
