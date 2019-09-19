function V=computeV(geo,angles,alphablocks,orig_index)

V=zeros(geo.nVoxel(1),geo.nVoxel(2) ,length(alphablocks),'single');
% V2=zeros(geo.nVoxel(1),geo.nVoxel(2) ,length(alphablocks),'single');
geo=checkGeo(geo,angles);

if ~isfield(geo,'mode')||~strcmp(geo.mode,'parallel')
    for ii=1:length(alphablocks)
        [x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
            -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
        auxang=alphablocks{ii};
        auxindex=orig_index{ii};
%         A = permute(auxang(1,:)+pi/2, [1 3 2]);
        for jj=1:length(auxang(1,:))
            V(:,:,ii)=V(:,:,ii)+single(((geo.DSO(auxindex(jj)) ./ (geo.DSO(auxindex(jj)) +  y.*sin(-auxang(1,jj)-pi/2) -  x.*cos(-auxang(1,jj)-pi/2))).^2).');
        end
%         V(:,:,ii)= sum(permute(single((geo.DSO ./ (geo.DSO +  y.*sin(-A) -  x.*cos(-A))).^2),[2 1 3]),3);
    end
else
    for ii=1:length(alphablocks)
        V(:,:,ii)=ones(geo.nVoxel(1:2).','single')*length(alphablocks{ii});
    end
end
end