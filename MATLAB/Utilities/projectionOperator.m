function op=projectionOperator(geo,acquisition,varargin)
%PROJECTIONOPERATOR Bind conventional angles or an explicit detector surface.
% op.forward(volume), op.adjoint(data), op.subset(oneBasedViewIndices).
% MATLAB volumes are xyz; surface data are rows x columns x views.
% Conventional Ax also uses rows(v) x columns(u) x views. Device is zero based.
p=inputParser; addParameter(p,'Device',0); parse(p,varargin{:}); device=p.Results.Device;
validateattributes(device,{'numeric'},{'scalar','integer','nonnegative','finite'});
op.volume_shape=double(geo.nVoxel(:))'; op.device=device;
validateattributes(op.volume_shape,{'double'},{'numel',3,'integer','positive','finite'});
if isstruct(acquisition)
    surface=acquisition; n=size(surface.source_positions_mm,2);
    assert(n>0,'TIGRE:operator:views','Acquisition must have views');
    op.projection_shape=[size(surface.pixel_centers_mm,1) size(surface.pixel_centers_mm,2) n];
    mask=true(op.projection_shape);
    if isfield(surface,'valid_mask')
        mask=surface.valid_mask;
        assert(islogical(mask),'TIGRE:operator:mask','Mask must be logical');
        if size(mask,3)==1, mask=repmat(mask,1,1,n); end
        assert(isequal(size3(mask),op.projection_shape),'TIGRE:operator:mask','Invalid mask shape');
    end
    op.valid_mask=mask;
    op.forward=@(x) Ax_surface(checked(x,op.volume_shape,[]),geo,surface,device);
    op.adjoint=@(b) Atb_surface(checked(b,op.projection_shape,mask),geo,surface,device);
    op.subset=@(idx) surfaceSubset(geo,surface,idx,n,device);
else
    angles=acquisition; if isvector(angles), angles=reshape(angles,1,[]); end
    geo=checkGeo(geo,angles); n=size(angles,2);
    ids=GpuIds(); ids.devices=int32(device);
    op.projection_shape=[double(geo.nDetector(2)) double(geo.nDetector(1)) n]; op.valid_mask=[];
    op.forward=@(x) Ax(checked(x,op.volume_shape,[]),geo,angles,'Siddon','gpuids',ids);
    op.adjoint=@(b) Atb(checked(b,op.projection_shape,[]),geo,angles,'matched','gpuids',ids);
    op.subset=@(idx) conventionalSubset(geo,angles,idx,n,device);
end
end
function out=surfaceSubset(geo,s,idx,n,device)
indices(idx,n); s.source_positions_mm=s.source_positions_mm(:,idx);
if size(s.pixel_centers_mm,4)>1, s.pixel_centers_mm=s.pixel_centers_mm(:,:,:,idx); end
for field={'valid_mask','detector_valid_mask'}
    name=field{1}; if isfield(s,name) && size(s.(name),3)>1, s.(name)=s.(name)(:,:,idx); end
end
% Beam effects are already represented by valid_mask; reconstruction uses
% exact rays rather than editable beam recipes.
if isfield(s,'source_beam'), s=rmfield(s,'source_beam'); end
out=projectionOperator(geo,s,'Device',device);
end
function out=conventionalSubset(geo,angles,idx,n,device)
indices(idx,n);
for field={'DSO','DSD','COR','offOrigin','offDetector','rotDetector'}
    name=field{1}; if isfield(geo,name) && size(geo.(name),2)==n, geo.(name)=geo.(name)(:,idx); end
end
out=projectionOperator(geo,angles(:,idx),'Device',device);
end
function indices(idx,n)
validateattributes(idx,{'numeric'},{'vector','nonempty','integer','>=',1,'<=',n,'finite'});
end
function x=checked(x,shape,mask)
assert(isnumeric(x)&&isreal(x)&&isequal(size3(x),shape),'TIGRE:operator:shape','Invalid array shape/type');
x=single(x); if ~isempty(mask), x(~mask)=0; end
assert(all(isfinite(x(:))),'TIGRE:operator:finite','Nonfinite valid samples');
end
function s=size3(x)
s=[size(x,1) size(x,2) size(x,3)]; if ndims(x)>3, s=[]; end
end
