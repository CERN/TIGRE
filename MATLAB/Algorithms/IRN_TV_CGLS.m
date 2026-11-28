function [x_outer,resL2,qualMeasOut]= IRN_TV_CGLS(proj,geo,angles,niter,varargin)
% IRN_TV_CGLS solves the CBCT problem with total variation regularization,
% by iteratively reweighted norms and CGLS (an inner-outer scheme).
%
%  IRN_TV_CGLS(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
%  IRN_TV_CGLS(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%  'lambda'  Value of regularization parameter lambda, default 10.
%  'Init'    Describes different initialization techniques.
%             * 'none'     : Initializes the image to zeros (default)
%             * 'FDK'      : initializes image to FDK reconstruction
%             * 'multigrid': Initializes image by solving the problem in
%                            small scale and increasing it when relative
%                            convergence is reached.
%             * 'image'    : Initialization using a user specified
%                            image. Not recommended unless you really
%                            know what you are doing.
%  'InitImg'    an image for the 'image' initialization. Avoid.
%  'groundTruth'  an image as ground truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
%  'restart'  true or false. By default the algorithm will restart when
%             loss of orthogonality is found.
%  'checkResidual' true or false. By default the algorithm computes the
%             residual ||Ax-b|| every iteration, which costs one extra
%             forward projection, and restarts or stops when it grows (loss
%             of orthogonality). false skips both, for speed; the residual
%             is still computed if it is asked for as an output.
%  'niter_outer'  Number of outer (reweighting) iterations, default 10.
%                 NITER is shared between them.
%  'outerRestart' What each outer iteration starts from:
%             * 'cold' : the initial image (default)
%             * 'warm' : the image of the previous outer iteration
%  'weights'      TV weights for the first outer iteration, of size
%                 geo.nVoxel', instead of computing them from the initial
%                 image (for example from a previous solve).
%
% OUTPUTS:
%
%    [img]                       will output the reconstructed image
%    [img,resL2]                 will also output the L2 norm of the
%                                residual, ||Ax-b||, after each iteration
%    [img,resL2,qualMeas]        will also output the quality measurements
%                                asked by the input 'QualMeas'
%
%  See also IRN_PICCS_CGLS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD.
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Malena Sabate Landman, Ander Biguri
%--------------------------------------------------------------------------

% IRN_PICCS_CGLS without a prior term is this algorithm: one implementation for both.
out=cell(1,max(nargout,1));
[out{:}]=IRN_PICCS_CGLS(proj,geo,angles,niter,[],varargin{:});
x_outer=out{1};
if nargout>1
    resL2=out{2};
end
if nargout>2
    qualMeasOut=out{3};
end
end
