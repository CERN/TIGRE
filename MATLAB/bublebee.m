
InitTIGRE;
folder='D:\CTdata\bumblebee_projections\20130430_HMX_446_RB_4fpp_40kVp_Mo_Bumblebee';

[proj,geo,angles]=NikonDataLoader(folder,'num_angles',100);


geo.nVoxel(3)=100;
geo.sVoxel(3)=geo.nVoxel(3).*geo.dVoxel(3);
geo.COR=17.1848*geo.dDetector(1);
geo.COR= 1.03014;

img=FDK(proj,geo,angles);
%%
tosave=squeeze(img(round(geo.nVoxel(1)/2),:,:));
tosave=max(tosave,0)./0.1;
% plotImg(img,'dim',3,'slice',50,'clims',[0, 0.1])
imwrite(tosave,"test.png")