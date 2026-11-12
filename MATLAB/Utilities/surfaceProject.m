function result=surfaceProject(back,data,geo,surface,device)
% Internal layout conversion for the shared native CUDA surface kernels.
if nargin<5, device=0; end
validateattributes(device,{'numeric'},{'real','scalar','integer','nonnegative','finite','<=',double(intmax('int32'))});
validateattributes(data,{'numeric'},{'real','nonsparse'});
assert(isstruct(geo) && all(isfield(geo,{'nVoxel','dVoxel'})), ...
    'TIGRE:surface:geometry','Volume geometry needs nVoxel and dVoxel.');
n=double(geo.nVoxel(:)); d=double(geo.dVoxel(:));
validateattributes(n,{'double'},{'numel',3,'real','integer','positive','finite','<=',double(intmax('int32'))});
validateattributes(d,{'double'},{'numel',3,'real','positive','finite'});
o=zeros(3,1); if isfield(geo,'offOrigin'), o=geo.offOrigin; end
validateattributes(o,{'numeric'},{'numel',3,'real','finite'});
assert(isstruct(surface) && all(isfield(surface,{'source_positions_mm','pixel_centers_mm'})), ...
    'TIGRE:surface:geometry','Surface needs source_positions_mm and pixel_centers_mm.');
s=surface.source_positions_mm; p=surface.pixel_centers_mm;
validateattributes(s,{'numeric'},{'2d','nrows',3,'nonempty','real','nonsparse','finite'});
validateattributes(p,{'numeric'},{'nonempty','real','nonsparse','finite'});
assert(size(p,3)==3 && ndims(p)<=4,'TIGRE:surface:pixels','Pixels must be rows x columns x 3 x views.');
v=size(s,2); r=size(p,1); c=size(p,2);
assert(size(p,4)==1 || size(p,4)==v,'TIGRE:surface:views','Pixel views must be one or match sources.');
m=true(r,c); if isfield(surface,'valid_mask'), m=surface.valid_mask; end
validateattributes(m,{'logical','numeric'},{'real','nonsparse','finite'});
assert(size(m,1)==r && size(m,2)==c && ...
    (size(m,3)==1 || size(m,3)==v) && ndims(m)<=3 && all(m(:)==0 | m(:)==1), ...
    'TIGRE:surface:mask','Mask must be binary rows x columns x (one or views).');
expected=n'; if back, expected=[r c v]; end
assert(size(data,1)==expected(1) && size(data,2)==expected(2) && ...
    size(data,3)==expected(3) && ndims(data)<=3,'TIGRE:surface:shape','Data shape does not match geometry.');
data=single(data);
if back
    valid=repmat(logical(m),1,1,v/size(m,3));
    data(~valid)=0; % Missing samples are absent, even if their stored values are NaN.
    data=permute(data,[2 1 3]);
end
assert(all(isfinite(data(:))),'TIGRE:surface:finite','Valid data must be finite float32.');
result=Surface_mex(single(back),data,n,single(d),single(o), ...
    single(s),single(permute(p,[3 2 1 4])),uint8(permute(m,[2 1 3])),single(device));
if ~back, result=permute(result,[2 1 3]); end
end
