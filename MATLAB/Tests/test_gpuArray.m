function test_gpuArray()
% Ax and Atb with gpuArray inputs (Parallel Computing Toolbox) give the same
% result as with host arrays, as gpuArrays on the same GPU. Needs a GPU, the
% toolbox, and TIGRE compiled with the toolbox installed.
if ~(exist('gpuArray','file') && gpuDeviceCount>0)
    fprintf('SKIP gpuArray: no Parallel Computing Toolbox or no GPU\n');
    return
end
for mode={'cone','parallel'}
    if strcmp(mode{1},'cone')
        nDetector=[80;60];
    else
        nDetector=[64;48]; % parallel beam: as nVoxel(2:3)
    end
    geo=defaultGeometry('nVoxel',[64;64;48],'nDetector',nDetector,'mode',mode{1});
    angles=linspace(0,2*pi,37);
    angles=angles(1:end-1);
    img=single(headPhantom(geo.nVoxel));
    for ptype={'Siddon','interpolated'}
        p=Ax(img,geo,angles,ptype{1});
        pg=Ax(gpuArray(img),geo,angles,ptype{1});
        check(pg,p,sprintf('Ax %s %s',mode{1},ptype{1}));
    end
    backtypes={'FDK','matched'};
    if strcmp(mode{1},'parallel'), backtypes={'FDK'}; end
    for btype=backtypes
        b=Atb(p,geo,angles,btype{1});
        bg=Atb(gpuArray(p),geo,angles,btype{1});
        check(bg,b,sprintf('Atb %s %s',mode{1},btype{1}));
    end
    % one projection: 2D projection data
    check(Ax(gpuArray(img),geo,angles(1)),Ax(img,geo,angles(1)),sprintf('Ax %s one angle',mode{1}));
    check(Atb(gpuArray(p(:,:,1)),geo,angles(1)),Atb(p(:,:,1),geo,angles(1)),sprintf('Atb %s one angle',mode{1}));
end
% 'gpuids' naming another GPU: host copy on that GPU, still a gpuArray out (as in Python)
if gpuDeviceCount>1
    g=gpuDevice();
    other=GpuIds();
    other.devices=int32(mod(g.Index,gpuDeviceCount)); % the next GPU, 0-based
    check(Ax(gpuArray(img),geo,angles,'gpuids',other),Ax(img,geo,angles),'Ax gpuids on another GPU');
end
% errors stay clear: wrong type
try
    Ax(gpuArray(double(img)),geo,angles);
    error('test_gpuArray:noError','double gpuArray should have been rejected');
catch err
    assert(strcmp(err.identifier,'TIGRE:Ax:InvalidInput'),'double gpuArray: wrong error %s',err.identifier);
end
fprintf('PASS gpuArray: Ax/Atb on gpuArrays match host arrays\n');
end

function check(got,expected,name)
assert(isa(got,'gpuArray'),'%s: result should be a gpuArray',name);
assert(strcmp(classUnderlying(got),'single'),'%s: result should be single',name);
assert(isequal(size(got),size(expected)),'%s: size %s, expected %s',name,mat2str(size(got)),mat2str(size(expected)));
got=gather(got);
err=norm(got(:)-expected(:))/max(norm(expected(:)),eps('single'));
assert(err<1e-5,'%s: relative difference %g',name,err);
fprintf('  %s: relative difference %.2g\n',name,err);
end
