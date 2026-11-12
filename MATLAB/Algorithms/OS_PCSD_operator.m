function [image,residual,quality]=OS_PCSD_operator(data,op,niter,varargin)
% Explicit opt-in OS_PCSD; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('OS_PCSD',data,op,niter,varargin{:});
end
