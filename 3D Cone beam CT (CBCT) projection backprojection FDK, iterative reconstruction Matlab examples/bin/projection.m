function proj2d = projection(data3d,param, iview)

angle_rad = param.deg(iview)/360*2*pi;
proj2d = (zeros(param.nu,param.nv,'single'));

[uu,vv] = meshgrid(param.us,param.vs);
[xx,yy] = meshgrid(param.xs,param.ys);

if param.gpu == 1
    data3d = gpuArray(single(data3d));
    rx = gpuArray(((xx.*cos(angle_rad) - yy.*sin(angle_rad)) - xx(1,1))/param.dx + 1);
    ry = gpuArray(((xx.*sin(angle_rad) + yy.*cos(angle_rad)) - yy(1,1))/param.dy + 1);
else
    rx = (((xx.*cos(angle_rad) - yy.*sin(angle_rad)) - xx(1,1))/param.dx + 1);
    ry = (((xx.*sin(angle_rad) + yy.*cos(angle_rad)) - yy(1,1))/param.dy + 1);
end

for iz = 1:param.nz   
    
    data3d(:,:,iz) = interp2(data3d(:,:,iz),rx,ry, param.interptype);
    
end

data3d(isnan(data3d))=0;

data3d = permute(data3d,[1 3 2]);

[xx,zz] = meshgrid(param.xs,param.zs);


for iy = 1:param.ny
    
    Ratio = (param.ys(iy)+param.DSO)/(param.DSD);
    
    pu = uu*Ratio;
    pv = vv*Ratio;    
    
    pu = (pu - xx(1,1))/(param.dx)+1; 
    pv = (pv - zz(1,1))/(param.dz)+1; 
    
    if param.gpu == 1
        tmp = gather(interp2(gpuArray(single(data3d(:,:,iy))),gpuArray(single(pv)),gpuArray(single(pu)),param.interptype));
    else
        tmp = (interp2((single(data3d(:,:,iy))),(single(pv)),(single(pu)),param.interptype));
    end
    
    tmp(isnan(tmp))=0;
    
    proj2d = proj2d + tmp';
end

dist = sqrt((param.DSD)^2 + uu.^2 + vv.^2)./(param.DSD)*param.dy;

proj2d = proj2d .* dist';





