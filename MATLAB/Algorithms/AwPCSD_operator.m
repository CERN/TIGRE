function [image,residual,quality]=AwPCSD_operator(data,op,niter,varargin)
% Explicit opt-in AwPCSD; conventional calls retain upstream behavior.
[image,residual,quality]=allOperatorRecon('AwPCSD',data,op,niter,varargin{:});
end
