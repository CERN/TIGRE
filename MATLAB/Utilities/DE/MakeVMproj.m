function [VM_proj] = MakeVMproj(material1_proj, material2_proj, energy, varargin)
%MAKEVMPROJ Synthesizes VM projections from basis material thickness
%projections.
%   
%   Coded by: Andrew Keeler
%   Contact: akeeler@luc.edu
%   Date: February 28, 2024

[material1, material2] = parse_inputs(varargin);

% read NIST database (values in TIGRE)
fid = fopen('./../../../Common/data/NIST_material_attenuation.json');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
NIST = jsondecode(str);

materials = fieldnames(NIST);

if ~isfield(NIST,material1)
    disp([material1, ' not a (yet) supported material in NIST dataset'])
    disp("Supported materials:")
    for i=1:len(materials)
        disp(materials{i})
    end
    error("")
end
if ~isfield(NIST,material2)
    disp([material2, ' not a (yet) supported material in NIST dataset'])
    disp("Supported materials:")
    for i=1:len(materials)
        disp(materials{i})
    end
    error("")
end

if (energy < max(min(NIST.(material1).energy),min(NIST.(material2).energy)) || energy > min(max(NIST.(material1).energy),max(NIST.(material2).energy)))
    warning('Requested energy is outside of expected range.  Results may be unreliable.')
end

%  Compute attenuation coefficients for basis materials at selected energy
material1_att = spline(NIST.(material1).energy, NIST.(material1).mu, energy);
material2_att = spline(NIST.(material2).energy, NIST.(material2).mu, energy);

VM_proj = material1_att*material1_proj + material2_att*material2_proj;
end



function [material1, material2]=parse_inputs(argin)

opts =  {'material1','material2'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:MakeVMproj:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option is not default, then extract value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end

        if isempty(ind)
            error('CBCT:MakeVMproj:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
        end
        val=argin{jj};
    end
    
    switch opt
        case 'material1'
            if default
                warning('Assuming Al for material1')
                material1='Al';
            else
                if ~ischar( val)
                    error('TIGRE:DEDecompose:InvalidInput','material1 must be a string')
                end
                material1=val;
             end
        case 'material2'
            if default
                warning('Assuming PMMA for material2')
                material2='PMMA';
            else
                if ~ischar( val)
                    error('TIGRE:DEDecompose:InvalidInput','material2 must be a string')
                end
                material2=val;
             end
              
        otherwise
            error('TIGRE:MakeVMproj:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in MakeVMproj()']);
    end
end
end
