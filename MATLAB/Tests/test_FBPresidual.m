function test_FBPresidual()
% FBP's second output is norm(proj - Ax(img)), also with a non-square
% detector. It used to subtract arrays in different layouts. Needs a GPU.
geo.DSD=200;
geo.DSO=100;
geo.nVoxel=[64;64;48];
geo.sVoxel=[64;64;48];
geo.dVoxel=geo.sVoxel./geo.nVoxel;
geo.nDetector=[64;48];
geo.dDetector=[1;1];
geo.sDetector=geo.nDetector.*geo.dDetector;
geo.offOrigin=[0;0;0];
geo.offDetector=[0;0];
geo.accuracy=0.5;
geo.mode='parallel';
angles=linspace(0,pi,60);
img=zeros(geo.nVoxel','single');
img(20:44,20:44,10:38)=1;
proj=Ax(img,geo,angles);
[res,errorL2]=FBP(proj,geo,angles);
residual=proj-Ax(res,geo,angles);
expected=norm(residual(:));
assert(abs(errorL2-expected)<=1e-4*expected,'FBP errorL2 should be norm(proj - Ax(res))');
fprintf('PASS FBP: errorL2 is the residual norm with a non-square detector\n');
end
