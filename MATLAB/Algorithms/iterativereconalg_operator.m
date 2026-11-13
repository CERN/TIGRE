function [image,residual,quality]=iterativereconalg_operator(data,op,niter,varargin)
% Explicit opt-in iterativereconalg; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('iterativereconalg',data,op,niter,varargin{:});
end
