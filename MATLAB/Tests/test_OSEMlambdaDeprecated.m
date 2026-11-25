function test_OSEMlambdaDeprecated()
% OSEM's 'lambda' and 'lambda_red' never changed the result. They still
% work, with a TIGRE:deprecated warning. Needs a GPU.
geo=defaultGeometry('nVoxel',[32;32;32]);
angles=linspace(0,2*pi,30);
proj=Ax(headPhantom(geo.nVoxel),geo,angles);
for opt={'lambda','lambda_red'}
    lastwarn('');
    OSEM(proj,geo,angles,1,opt{1},0.5,'Verbose',0);
    [~,id]=lastwarn();
    assert(strcmp(id,'TIGRE:deprecated'),['OSEM should warn that ''' opt{1} ''' is deprecated']);
end
lastwarn('');
OSEM(proj,geo,angles,1,'Verbose',0);
[~,id]=lastwarn();
assert(~strcmp(id,'TIGRE:deprecated'),'OSEM should not warn without the deprecated options');
fprintf('PASS OSEM: lambda and lambda_red give a deprecation warning\n');
end
