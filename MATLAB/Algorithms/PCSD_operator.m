function [image,residual,quality]=PCSD_operator(data,op,niter,varargin)
% Explicit opt-in PCSD; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('PCSD',data,op,niter,varargin{:});
end
