function [image,residual,quality]=LSMR_operator(data,op,niter,varargin)
% Explicit opt-in LSMR; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('LSMR',data,op,niter,varargin{:});
end
