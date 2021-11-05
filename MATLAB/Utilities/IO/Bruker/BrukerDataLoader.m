function [proj,geo,angles]=BrukerDataLoader(filepath,varargin)
%BrukerDataLoader(filepath) Loads Bruker uCT datasets into TIGRE standard
%
%   BrukerDataLoader(filepath OPT,VAL,...) uses options and values.
%      These are options in case you don't want to load the entire
%      dataset, but only particular sets of projections.
%      The possible options in OPT are:
%           'sampling': type of sampling. default 'equidistant' Can be:
%                  'equidistant': equidistantly sample the entire set
%                       of angles. 'num_angles' should be set for partial
%                       loading of the data.
%                  'step': sample the entire set of projections every
%                         'sampling_step' angles.
%                  'continous': Load the first 'num_angles' amount of
%                             angles only.
%
%           'num_angles': Number of total angles to load. Default all of
%                    them. Useful for 'equidistant' and 'continous' loading
%
%           'sampling_step': step to load when loading projections.
%                   Default=1. Useful for 'step' loading.
%
%           'dataset_number': some Bruker data has several datasets in a
%                             single folder. This can be the dataset number
%                             of 'all' to load all of them
%% Parse inputs

dataset_number=-1;
for ii=1:2:length(varargin)
    if strcmpi('dataset_number',varargin{ii})
        dataset_number=varargin{ii+1};
        break;
    end
end

if strcmpi(dataset_number,'all')
    num_scans=find_number_of_datasets(filepath);
    if ~isempty(num_scans)
        disp("Loading all scans in folder, assuming the same geometry for all")
        
        angles=[];
        proj=[];
        for ii=1:num_scans
             disp(['Loading scan number ' , num2str(ii) , '...'])
             [geo,auxangles]=readBrukerGeometry(filepath,ii-1);
             [auxproj,geo,auxangles]=loadBrukerProjections(filepath,geo,auxangles,varargin{:});
             angles=[angles auxangles];
             proj=cat(3,proj,auxproj);
             disp('');
        end
        return;
    else
        [geo,angles]=readBrukerGeometry(filepath);
        [proj,geo,angles]=loadBrukerProjections(filepath,geo,angles,varargin{:});
        return;
    end
else
    [geo,angles]=readBrukerGeometry(filepath,dataset_number);
    [proj,geo,angles]=loadBrukerProjections(filepath,geo,angles,varargin{:});
    return;
end
end

function num=find_number_of_datasets(folder_or_file)
file = dir([folder_or_file,'/*.log']); %
if isempty(file)
    error(['No .log file found in folder: ', folder_or_file]);
end

filename=file(1).name;
fid=fopen([folder_or_file,'/',filename]);
xtekctText = textscan(fid, '%s %s', 'Delimiter', '=', 'HeaderLines', 1, 'CommentStyle', '[');
fclose(fid);

num=str2double(xtekctText{2}(strcmp('Number of connected scans', xtekctText{1})));

end