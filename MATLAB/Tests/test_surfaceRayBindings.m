function report=test_surfaceRayBindings(outputFile)
% Independent small-matrix CUDA binding checks; needs a compiled Surface_mex.
% Does not require GUI, Python, extra toolboxes, geometry presets or solvers.
if nargin<1, outputFile=''; end
rng(7); cases=0;
geo=struct('nVoxel',[6;5;4],'dVoxel',[1.2;.8;1.5],'offOrigin',[.13;-.21;.17]);
x=reshape(single(0:119)/120,geo.nVoxel');
for kind=1:3
    surface=fixture(kind);
    A=rayMatrix(geo,surface);
    q=reshape(single(sin(1:numel(surface.valid_mask))),size(surface.valid_mask));
    ax=Ax_surface(x,geo,surface); at=Atb_surface(q,geo,surface);
    close(ax(:),A*double(x(:))); close(at(:),A'*double(q(:)));
    left=double(ax(:))'*double(q(:)); right=double(x(:))'*double(at(:));
    assert(abs(left-right)<2e-5*max([1 abs(left) abs(right)]));
    assert(all(ax(~surface.valid_mask)==0)); cases=cases+1;
end
% Same coordinates/data also form an optional Python cross-binding fixture.
if ~isempty(outputFile), save(outputFile,'geo','surface','x','q','ax','at','-v7'); end
dirty=q; dirty(~surface.valid_mask)=NaN; before=dirty;
close(Atb_surface(dirty,geo,surface),at); assert(isequaln(dirty,before)); cases=cases+1;
bad=q; bad(find(surface.valid_mask,1))=NaN;
expectError(@() Atb_surface(bad,geo,surface),'TIGRE:surface:finite'); cases=cases+1;

% Shared pixel table with per-view masks, then shared table and mask.
fixed=surface; fixed.pixel_centers_mm=surface.pixel_centers_mm(:,:,:,1);
fixed.valid_mask(:,:,2)=true; A=rayMatrix(geo,fixed);
close(Ax_surface(x,geo,fixed),reshape(A*double(x(:)),size(q)));
close(Atb_surface(q,geo,fixed),reshape(A'*double(q(:)),geo.nVoxel')); cases=cases+1;
fixed.valid_mask=fixed.valid_mask(:,:,1); A=rayMatrix(geo,fixed);
close(Ax_surface(x,geo,fixed),reshape(A*double(x(:)),size(q))); cases=cases+1;

for axis=1:3
    s=geo.offOrigin; p=s; s(axis)=30; p(axis)=-30;
    one=struct('source_positions_mm',s,'pixel_centers_mm',reshape(p,1,1,3));
    value=Ax_surface(ones(geo.nVoxel','single'),geo,one);
    assert(abs(double(value)-geo.nVoxel(axis)*geo.dVoxel(axis))<2e-5); cases=cases+1;
end
% Segment endpoints inside, diagonal ties, lower/upper faces, outside and zero ray.
box=struct('nVoxel',[4;4;4],'dVoxel',ones(3,1),'offOrigin',zeros(3,1));
starts=[5 0 0;0 0 0;.2 .3 .4;5 5 5;5 -2 0;5 2 0;5 3 0;0 0 0];
ends=[0 0 0;-5 0 0;.4 .5 .6;-5 -5 -5;-5 -2 0;-5 2 0;-5 3 0;0 0 0];
volume=reshape(single(0:63)/64,[4 4 4]);
for k=1:size(starts,1)
    one=struct('source_positions_mm',starts(k,:)','pixel_centers_mm',reshape(ends(k,:),1,1,3));
    A=rayMatrix(box,one);
    close(Ax_surface(volume,box,one),A*double(volume(:)));
    close(Atb_surface(single(1),box,one),reshape(A',[4 4 4])); cases=cases+1;
end

for k=1:3
    expectError(@() Ax_surface(x,geo,surface,99999),'TIGRE:surface:CUDA');
    close(Ax_surface(x,geo,surface),ax); cases=cases+1;
end
expectError(@() Ax_surface(x(:,:,1),geo,surface),'TIGRE:surface:shape'); cases=cases+1;
invalid=surface; invalid.valid_mask=2*ones(size(q));
expectError(@() Ax_surface(x,geo,invalid),'TIGRE:surface:mask'); cases=cases+1;
invalid=surface; invalid.pixel_centers_mm=complex(invalid.pixel_centers_mm,1);
expectError(@() Ax_surface(x,geo,invalid),''); cases=cases+1;
invalidGeo=geo; invalidGeo.dVoxel(1)=1e-100;
expectError(@() Ax_surface(x,invalidGeo,surface),'TIGRE:surface:geometry'); cases=cases+1;

% Malformed direct packed calls must raise before native memory/launch work.
args={single(0),x,double(geo.nVoxel),single(geo.dVoxel),single(geo.offOrigin), ...
    single(surface.source_positions_mm),single(permute(surface.pixel_centers_mm,[3 2 1 4])), ...
    uint8(permute(surface.valid_mask,[2 1 3])),single(0)};
badargs=args; badargs{7}=zeros(3,0,1,'single'); badargs{8}=zeros(0,1,'uint8');
expectError(@() Surface_mex(badargs{:}),'TIGRE:surface:views'); cases=cases+1;
badargs=args; badargs{3}=repmat(double(intmax('int32')),3,1);
expectError(@() Surface_mex(badargs{:}),'TIGRE:surface:size'); cases=cases+1;
badargs=args; badargs{8}(:)=2;
expectError(@() Surface_mex(badargs{:}),'TIGRE:surface:mask'); cases=cases+1;
badargs=args; badargs{2}=single(0);
expectError(@() Surface_mex(badargs{:}),'TIGRE:surface:data'); cases=cases+1;
badargs=args; badargs{6}(1)=Inf;
expectError(@() Surface_mex(badargs{:}),'TIGRE:surface:finite'); cases=cases+1;
close(Ax_surface(x,geo,surface),ax); cases=cases+1;
report=struct('cases',cases,'matlab',version,'mex',which('Surface_mex'));
fprintf('PASS MATLAB surface ray bindings: %d cases; MEX %s\n',cases,report.mex);
end

function surface=fixture(kind)
beta=[.1 .7 1.6]; pixels=zeros(3,5,3,3); source=zeros(3,3);
u=linspace(-7,7,5); local=zeros(3,5,3);
local(:,:,1)=repmat(-20+double(kind==2)*.03*u.^2,3,1);
local(:,:,2)=repmat(u,3,1); local(:,:,3)=repmat(linspace(-4,4,3)',1,5);
if kind==3
    local(:,:,1)=repmat([-7 -19 -20 -19 -7],3,1);
    local(:,:,2)=repmat([-4 -4 0 4 4],3,1);
end
for v=1:3
    c=cos(beta(v)); s=sin(beta(v)); R=[c -s 0;s c 0;0 0 1];
    source(:,v)=R*[30;0;.3];
    pixels(:,:,:,v)=reshape(reshape(local,[],3)*R',[3 5 3]);
end
mask=true(3,5,3); mask(:,3,:)=false;
surface=struct('source_positions_mm',single(source),'pixel_centers_mm',single(pixels),'valid_mask',mask);
end

function A=rayMatrix(geo,surface)
% Intersect each voxel box separately; do not reproduce Siddon's traversal.
n=double(geo.nVoxel(:)); spacing=double(geo.dVoxel(:));
origin=double(geo.offOrigin(:))-.5*n.*spacing;
rows=size(surface.pixel_centers_mm,1); cols=size(surface.pixel_centers_mm,2);
views=size(surface.source_positions_mm,2); A=zeros(rows*cols*views,prod(n));
mask=true(rows,cols,views);
if isfield(surface,'valid_mask'), mask=repmat(logical(surface.valid_mask),1,1,views/size(surface.valid_mask,3)); end
for sample=1:size(A,1)
    if ~mask(sample), continue; end
    [r,c,v]=ind2sub([rows cols views],sample);
    pv=min(v,size(surface.pixel_centers_mm,4));
    source=double(surface.source_positions_mm(:,v));
    ray=double(reshape(surface.pixel_centers_mm(r,c,:,pv),3,1))-source;
    lengthMm=norm(ray);
    for voxel=1:size(A,2)
        [i,j,k]=ind2sub(n',voxel); lo=origin+([i;j;k]-1).*spacing; hi=lo+spacing;
        enter=0; leave=1;
        for axis=1:3
            if ray(axis)==0
                if source(axis)<lo(axis) || source(axis)>=hi(axis), leave=-1; break; end
            else
                limits=sort(([lo(axis);hi(axis)]-source(axis))/ray(axis));
                enter=max(enter,limits(1)); leave=min(leave,limits(2));
            end
        end
        A(sample,voxel)=max(0,leave-enter)*lengthMm;
    end
end
end

function close(actual,expected)
actual=double(actual(:)); expected=double(expected(:));
assert(numel(actual)==numel(expected) && all(abs(actual-expected)<=2e-5+2e-4*abs(expected)), ...
    'TIGRE:test:reference','Disagreement with independent reference.');
end

function expectError(fn,identifier)
try
    result=fn(); %#ok<NASGU>
catch err
    assert(isempty(identifier) || strcmp(err.identifier,identifier), ...
        'TIGRE:test:error','Expected %s; got %s: %s',identifier,err.identifier,err.message);
    return
end
error('TIGRE:test:missingError','Expected an exception.');
end
