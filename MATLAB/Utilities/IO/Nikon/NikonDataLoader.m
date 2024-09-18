function [proj,geo,angles]=NikonDataLoader(filepath,varargin)
%NikonDataLoader(filepath) Loads Nikon uCT datasets into TIGRE standard
%
%   NikonDataLoader(filepath OPT,VAL,...) uses options and values. 
%      These are options in case you don't want to load the entire
%      dataset, but only particular sets of projections. 
%      The possible options in OPT are:
%           'sampling': type of sampling. default 'equidistant' Can be:
%                  'equidistant': equidistantly sample the entire set
%                       of angles. 'num_angles' should be set for partial
%                       loading of the data. 
%                  'step': sample the entire set of projections every
%                         'sampling_step' angles.
%                  'continuous': Load the first 'num_angles' amount of
%                             angles only.
%
%           'num_angles': Number of total angles to load. Default all of
%                    them. Useful for 'equidistant' and 'continuous' loading
%
%           'sampling_step': step to load when loading projections. 
%                   Default=1. Useful for 'step' loading.
%%

[geo,angles]=readXtekctGeometry(filepath);
[proj,geo,angles]=loadNikonProjections(filepath,geo,angles,varargin{:});


end