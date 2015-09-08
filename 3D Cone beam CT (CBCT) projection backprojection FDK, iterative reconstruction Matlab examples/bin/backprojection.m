function vol = backprojection(proj,param,iview)

angle_rad = param.deg(iview)/360*2*pi;
vol = zeros(param.nx,param.ny,param.nz,'single');

[xx,yy] = meshgrid(param.xs,param.ys);

rx = xx.*cos(angle_rad-pi/2) + yy.*sin(angle_rad-pi/2);
ry = -xx.*sin(angle_rad-pi/2) + yy.*cos(angle_rad-pi/2);

pu = single(((rx.*(param.DSD)./(ry + param.DSO))+param.us(1))/(-param.du) + 1);
Ratio = (single(param.DSO.^2./(param.DSO+ry).^2));    
if param.gpu == 1
    pu = gpuArray(single(pu));
    proj = gpuArray(single(proj));
    Ratio = gpuArray(Ratio);
end    

for iz = 1:param.nz   
    if param.gpu == 1
        pv = gpuArray(single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1));
        vol(:,:,iz) = gather(Ratio.*interp2(proj',pu,pv,param.interptype));
    else
        pv = single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1);
        vol(:,:,iz) = (Ratio.*interp2(proj',pu,pv,param.interptype));
    end 
end

vol(isnan(vol))=0;

return

















