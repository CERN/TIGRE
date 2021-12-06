function V=computeV(geo,angles,alphablocks,orig_index)

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
        
        V(:,:,ii) = mean(Atb(ones(geo.nDetector(2),geo.nDetector(1),length(auxang),'single'),auxgeo,auxang),3);
    end
    V(V==0.0)=Inf;
else
    for ii=1:length(alphablocks)
        V(:,:,ii)=ones(geo.nVoxel(1:2).','single')*length(alphablocks{ii});
    end
end
end
