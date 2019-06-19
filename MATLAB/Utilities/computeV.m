function V=computeV(geo,angles,alphablocks)

V=zeros(geo.nVoxel(1),geo.nVoxel(2) ,length(alphablocks),'single');
if ~isfield(geo,'mode')||~strcmp(geo.mode,'parallel')
    for ii=1:length(alphablocks)
        [x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
            -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
        auxang=alphablocks{ii};
        A = permute(auxang(1,:)+pi/2, [1 3 2]);
        V(:,:,ii)= sum(permute(single((geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2),[2 1 3]),3);
    end
else
    for ii=1:length(alphablocks)
        V(:,:,ii)=ones([geo.nVoxel(1:2).'],'single')*length(alphablocks{ii});
    end
end
end