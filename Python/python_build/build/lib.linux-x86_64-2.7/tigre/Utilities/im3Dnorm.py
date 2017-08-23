import numpy as np

def im3DNORM(img,normind,varargin=None):
    '''
    % IMAGE3DNORM computes the desired image norm
    %   IMAGE3DNORM(IMG,NORMIND) computes the norm if image IMG using the norm
    %   defined in NORMING
    %
    %   IMG         A 3D image
    %   NORMIND     {non-zero int, inf, -inf, 'fro', 'nuc'}
    %               'TV': TV norm
    %
    %
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
    %                     https://github.com/CERN/TIGRE/license.txt
    %
    % Contact:            tigre.toolbox@gmail.com
    % Codes:              https://github.com/CERN/TIGRE/
    % Coded by:           Ander Biguri
    %--------------------------------------------------------------------------
    '''

    if normind is [np.inf,-np.inf,'fro','nuc']:
        return np.linalg.norm(img.ravel(),normind)
    if type(normind) is int:
        return np.linalg.norm(img.ravel(),normind)
    if normind is 'TV':
        raise Exception('TV not yet implemented')


