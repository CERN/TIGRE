% Small finite-source detector example; no reconstruction or gap completion.
% Run Compile from MATLAB/ first, then add Utilities/ and the MEX output path.
geo=struct('nVoxel',[16;16;8],'dVoxel',ones(3,1),'offOrigin',zeros(3,1));
volume=ones(geo.nVoxel','single'); source=[40;0;0];
theta=linspace(-.3,.3,24);
pixels=zeros(8,24,3);
pixels(:,:,1)=repmat(40-80*cos(theta),8,1);
pixels(:,:,2)=repmat(80*sin(theta),8,1);
pixels(:,:,3)=repmat(((0:7)-3.5)',1,24);
curved=struct('source_positions_mm',source,'pixel_centers_mm',pixels);

local=zeros(8,12,3); local(:,:,2)=repmat(linspace(-6,6,12),8,1);
local(:,:,3)=repmat(((0:7)-3.5)',1,12);
pixels=zeros(8,38,3); mask=false(8,38);
angles=[-pi/2 0 pi/2]; centres=[-34 -6.5 0;-40 0 0;-34 6.5 0];
for panel=1:3
    c=cos(angles(panel)); s=sin(angles(panel)); R=[c -s 0;s c 0;0 0 1];
    columns=(1:12)+13*(panel-1);
    pixels(:,columns,:)=reshape(reshape(local,[],3)*R'+centres(panel,:),8,12,3);
    mask(:,columns)=true;
end
modular=struct('source_positions_mm',source,'pixel_centers_mm',pixels,'valid_mask',mask);
surfaces={curved,modular}; names={'curved','three-panel U'};
for example=1:numel(surfaces)
    projections=Ax_surface(volume,geo,surfaces{example},0);
    back=Atb_surface(projections,geo,surfaces{example},0);
    assert(all(isfinite(back(:))));
    fprintf('%s: [%s], projection norm %.6g\n',names{example},num2str(size(projections)),norm(projections(:)));
end
