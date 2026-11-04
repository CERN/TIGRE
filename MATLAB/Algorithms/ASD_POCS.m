function [ f,qualMeasOut ] = ASD_POCS(proj,geo,angles,maxiter,varargin)
%ASD_POCS Solves the ASD_POCS total variation constrained image in 3D
% tomography.
%
%   ASD_POCS(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
%   ASD_POCS(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%
%   'lambda':      Sets the value of the hyperparameter for the SART iterations.
%                  Default is 1
%
%   'lambdared':   Reduction of lambda. Every iteration
%                  lambda = lambdared * lambda. Default is 0.99
%
%   'init':        Describes different initialization techniques.
%                   •  'none'     : Initializes the image to zeros (default)
%
%                   •  'FDK'      : Initializes image to FDK reconstruction
%   'TViter':      Defines the amount of TV iterations performed per SART
%                  iteration. Default is 20
%
%   'alpha':       Defines the TV hyperparameter. Default is 0.002
%
%   'alpha_red':   Defines the reduction rate of the TV hyperparameter
%
%   'Ratio':       The maximum allowed image/TV update ratio. If the TV
%                  update changes the image more than this, the parameter
%                  will be reduced. Default is 0.95
%   'maxL2err'     Maximum L2 error to accept an image as valid. This
%                  parameter is crucial for the algorithm, determines at
%                  what point an image should not be updated further.
%                  Default is 20% of the FDK L2 norm.
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
%
% 'OrderStrategy'  Chooses the subset ordering strategy. Options are
%                  'ordered' : uses them in the input order, but divided
%                  'random'  : orders them randomly
%                  'angularDistance': chooses the next subset with the
%                                     biggest angular distance with the ones used.
% 'redundancy_weighting': true or false. Default is true. Applies data
%                         redundancy weighting to projections in the update step
%                         (relevant for offset detector geometry)
%  'groundTruth'  an image as ground truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
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
% Coded by:           Ander Biguri and Manasavee Lohvithee
%--------------------------------------------------------------------------



%% parse inputs
    if nargin < 4
        error('TIGRE:ASD_POCS:InvalidInput','Not enough inputs')
    end
    [ f,qualMeasOut ] = OS_ASD_POCS(proj,geo,angles,maxiter, varargin{:}, 'BlockSize', 1);
end
