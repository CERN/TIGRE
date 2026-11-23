function test_vectorGeometry()
% Conversion between TIGRE cone-beam geometries and per-view [source, detector centre, u, v] vectors.
% Needs no GPU. The fixed case at the end is shared with the Python tests, so both languages use one convention.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Utilities'));
tol=1e-9;

%% Circular scan: known vectors, and back to a rotation about z only
geo=baseGeo();
geo.DSO=300; geo.DSD=500;
a=linspace(0,2*pi,13); a=a(1:12);
vec=from_TIGRE_to_vector(geo,a);
z=zeros(size(a));
expected=[300*cos(a); 300*sin(a); z; -200*cos(a); -200*sin(a); z; -1.1*sin(a); 1.1*cos(a); z; z; z; z+1.3].';
assertClose(vec,expected,1e-12,'circular scan vectors');
[geo2,angles2]=from_vector_to_TIGRE(geo,vec);
assertClose(unwrap(angles2(1,:)),a,1e-12,'circular scan alpha');
assertClose(angles2(2:3,:),zeros(2,12),1e-12,'circular scan theta, psi');
assertClose(geo2.rotDetector,zeros(3,12),1e-12,'circular scan rotDetector');
assertClose(geo2.offDetector,zeros(2,12),1e-12,'circular scan offDetector');
assertClose(geo2.DSD,500*ones(1,12),1e-12,'circular scan DSD');

%% Random geometries: TIGRE -> vectors -> TIGRE describes the same scan
rng(0);
for trial=1:20
    n=7;
    geo=baseGeo();
    geo.DSO=250+100*rand(1,n); geo.DSD=geo.DSO+150+100*rand(1,n);
    geo.offOrigin=8*rand(3,n)-4; geo.offDetector=16*rand(2,n)-8;
    geo.rotDetector=2*rand(3,n)-1; geo.COR=6*rand(1,n)-3;
    angles=2*pi*rand(3,n)-pi;
    vec=from_TIGRE_to_vector(geo,angles);
    [geo2,angles2]=from_vector_to_TIGRE(geo,vec);
    assertClose(from_TIGRE_to_vector(geo2,angles2),vec,tol,'TIGRE round trip');
    assertClose(geo2.dDetector,geo.dDetector,1e-12,'dDetector');
    assert(isequal(size(geo2.offOrigin),[3 n]) && ~any(geo2.offOrigin(:)),'offOrigin should be zero');

    first=from_TIGRE_to_vector(geo,angles,'detector_reference','first_pixel');
    shift=vec(:,7:9)*(geo.nDetector(1)-1)/2+vec(:,10:12)*(geo.nDetector(2)-1)/2;
    assertClose(first(:,4:6),vec(:,4:6)-shift,tol,'first pixel');
    [geo3,angles3]=from_vector_to_TIGRE(geo,first,'detector_reference','first_pixel');
    assertClose(from_TIGRE_to_vector(geo3,angles3),vec,tol,'first pixel round trip');
end

%% Static linear source array and one fixed detector (no common rotation centre)
n=9;
src=[300*ones(n,1) linspace(-120,120,n).' linspace(-10,10,n).'];
vec=[src repmat([-200 5 -3 0 1.1 0 0 0 -1.3],n,1)];
[geo2,angles2]=from_vector_to_TIGRE(baseGeo(),vec);
assertClose(from_TIGRE_to_vector(geo2,angles2),vec,tol,'static source array');

%% Refusals
geo=baseGeo();
good=[300 0 0 -200 0 0 0 1.1 0 0 0 -1.3];
expectError(@() from_vector_to_TIGRE(geo,[good; good(1:9) 0 1.3*sin(0.01) -1.3*cos(0.01)]),'not perpendicular');
expectError(@() from_vector_to_TIGRE(geo,[good; good.*[1 1 1 1 1 1 1 1.001 1 1 1 1]]),'pixel size');
expectError(@() from_vector_to_TIGRE(geo,[0 0 0 good(4:12)]),'centre of the image');
expectError(@() from_vector_to_TIGRE(geo,[-200 50 0 good(4:12)]),'plane of the detector');
expectError(@() from_vector_to_TIGRE(geo,zeros(2,11)),'Nx12');
expectError(@() from_TIGRE_to_vector(geo,0,'detector_reference','corner'),'detector_reference');

%% Shared with Python (test_vector_geometry.py::test_matches_matlab_fixture)
geo=baseGeo();
geo.DSO=300; geo.DSD=480; geo.offOrigin=[1;-2;3]; geo.offDetector=[-6;4];
geo.rotDetector=[0.1;-0.2;0.3]; geo.COR=1.5;
vec=from_TIGRE_to_vector(geo,[0.7 2.5; 0.2 -0.1; -0.3 0.4]);
assertClose(vec,fixture(),1e-9,'fixture shared with Python');

fprintf('PASS vector geometry: round trips, known cases, refusals and the Python fixture\n');
end

function geo=baseGeo()
geo.mode='cone';
geo.nVoxel=[32;32;32]; geo.sVoxel=[64;64;64]; geo.dVoxel=[2;2;2];
geo.nDetector=[48;40]; geo.dDetector=[1.1;1.3]; geo.sDetector=geo.nDetector.*geo.dDetector;
end

function v=fixture()
v=[ 2.699824759324622e+02  1.162923904853032e+02 -5.993881829360623e+01 -1.621642623444666e+02 -6.973146434746344e+01  3.543582209769166e+01 -7.050992103620070e-01  8.371232783412641e-01  1.098167583115110e-01  8.360878861563337e-02 -9.926385137791681e-02  1.293505414861434e+00
   -2.920792266289853e+02  7.174692812602419e+01  2.458579979142952e+01  1.739759740279454e+02 -3.619112849413428e+01 -1.533820143203991e+01  9.139008720326046e-02 -1.095637089932312e+00  3.503168745057018e-02  3.499071258227064e-01  6.915757727574623e-02  1.250112888024055e+00];
end

function assertClose(actual,expected,tol,what)
assert(isequal(size(actual),size(expected)),'%s: size %s, expected %s',what,mat2str(size(actual)),mat2str(size(expected)));
err=max(abs(actual(:)-expected(:)));
assert(err<=tol,'%s: max error %g > %g',what,err,tol);
end

function expectError(f,fragment)
try
    f();
catch err
    assert(~isempty(strfind(err.message,fragment)),'Expected an error containing "%s", got "%s"',fragment,err.message);
    return
end
error('Expected an error containing "%s"',fragment);
end
