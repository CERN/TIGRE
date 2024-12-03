function W=computeW(geo,angles,gpuids)

    geoaux=geo;
    geoaux.sVoxel([1 2])=geo.sDetector(1)*1.1; % a bit bigger, to avoid numerical division by zero (small number)
    geoaux.sVoxel(3)=max(geo.sDetector(2),geo.sVoxel(3)); % make sure lines are not cropped. One is for when image is bigger than detector and viceversa
    geoaux.nVoxel=[2,2,2]'; % accurate enough?
    geoaux.dVoxel=geoaux.sVoxel./geoaux.nVoxel;
    W=Ax(ones(geoaux.nVoxel','single'),geoaux,angles,'Siddon','gpuids',gpuids);
    W(W<min(geo.dVoxel)/4)=Inf;
    W=1./W;
    W(W>0.1)=0.1;