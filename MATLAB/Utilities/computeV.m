function V=computeV(geo,angles,alphablocks,orig_index,varargin)
if nargin>=6
    [gpuids] = parse_inputs(varargin{1:length(varargin)});
elseif nargin==4
    gpuids = GpuIds();
else
    error('Wrong amount of inputs, 4 or 6 expected');
end
    
V=zeros(geo.nVoxel(1),geo.nVoxel(2) ,length(alphablocks),'single');
geo=checkGeo(geo,angles);

if ~isfield(geo,'mode')||~strcmp(geo.mode,'parallel')
    for ii=1:length(alphablocks)
        auxang=alphablocks{ii};
        auxindex=orig_index{ii};
        auxgeo = geo;
        % shrink the volume to avoiding zeros in backprojection
        auxgeo.sVoxel = auxgeo.sVoxel * max(auxgeo.sVoxel(1:2)/norm(auxgeo.sVoxel(1:2),2)) * 0.9;
        auxgeo.dVoxel = auxgeo.sVoxel ./ auxgeo.nVoxel;
        % subset of projection angles
        auxgeo.DSD = geo.DSD(auxindex);
        auxgeo.DSO = geo.DSO(auxindex);
        auxgeo.offOrigin = geo.offOrigin(:,auxindex);
        auxgeo.offDetector = geo.offDetector(:,auxindex);
        auxgeo.rotDetector = geo.rotDetector(:,auxindex);
        auxgeo.COR = geo.COR(auxindex);
                V(:,:,ii) = mean(Atb(ones(geo.nDetector(2),geo.nDetector(1),size(auxang,2),'single'),auxgeo,auxang, 'gpuids', gpuids),3);

        %V(:,:,ii) = mean(Atb(ones(geo.nDetector(2),geo.nDetector(1),length(auxang),'single'),auxgeo,auxang, 'gpuids', gpuids),3);
    end
else
    for ii=1:length(alphablocks)
        V(:,:,ii)=ones(geo.nVoxel(1:2).','single')*length(alphablocks{ii});
    end
end
end

function [gpuids]=parse_inputs(varargin)
    %fprintf('parse_inputs0(varargin (%d))\n', length(varargin));
    if isempty(varargin)
        gpuids = GpuIds();
    else
        % create input parser
        p=inputParser;
        % add optional parameters
        addParameter(p,'gpuids', GpuIds());
        %execute
        parse(p,varargin{:});
        %extract
        gpuids=p.Results.gpuids;
    end
end

