function V=computeV(geo,angles,alphablocks,orig_index, gpuids)

V=zeros(geo.nVoxel(1),geo.nVoxel(2) ,length(alphablocks),'single');
geo=checkGeo(geo,angles);

if ~isfield(geo,'mode')||~strcmp(geo.mode,'parallel')
    for ii=1:length(alphablocks)
        auxang=alphablocks{ii};
        auxindex=orig_index{ii};
        auxgeo = geo;
        % expand the detector to avoiding zeros in backprojection
%         maxsize=max(auxgeo.sVoxel+geo.offOrigin(:,auxindex),[],2);
%         auxgeo.sDetector=max(auxgeo.sDetector , [maxsize(1); maxsize(3)] *geo.DSD/geo.DSO);
%         auxgeo.dDetector = auxgeo.sDetector ./ auxgeo.nDetector;
        % subset of projection angles
        auxgeo.DSD = geo.DSD(auxindex);
        auxgeo.DSO = geo.DSO(auxindex);
        auxgeo.offOrigin = geo.offOrigin(:,auxindex);
        auxgeo.offDetector = geo.offDetector(:,auxindex);
        auxgeo.rotDetector = geo.rotDetector(:,auxindex);
        auxgeo.COR = geo.COR(auxindex);
        %auxgeo=geo;
        V(:,:,ii) = mean(Atb(ones(geo.nDetector(2),geo.nDetector(1),size(auxang,2),'single'),auxgeo,auxang, 'gpuids', gpuids),3)+0.000001;
    end
else
    for ii=1:length(alphablocks)
        V(:,:,ii)=ones(geo.nVoxel(1:2).','single')*length(alphablocks{ii});
    end
end
end
