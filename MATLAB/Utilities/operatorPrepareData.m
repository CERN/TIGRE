function b=operatorPrepareData(data,op)
%OPERATORPREPAREDATA Validate the complete contract before any fast return.
assert(isstruct(op),'TIGRE:operator:contract','Operator must be a struct');
for key={'volume_shape','projection_shape'}
    name=key{1}; assert(isfield(op,name),'TIGRE:operator:contract','Missing %s',name);
    validateattributes(op.(name),{'numeric'},{'vector','numel',3,'integer','positive','finite'});
end
for key={'forward','adjoint','subset'}
    name=key{1}; assert(isfield(op,name)&&isa(op.(name),'function_handle'), ...
        'TIGRE:operator:contract','Missing callable %s',name);
end
assert(isfield(op,'valid_mask'),'TIGRE:operator:contract','Missing valid_mask');
mask=op.valid_mask;
if ~isempty(mask)
    assert(islogical(mask)&&ndims(mask)<=3&& ...
        isequal([size(mask,1) size(mask,2) size(mask,3)],double(op.projection_shape(:))'), ...
        'TIGRE:operator:mask','Invalid boolean measurement mask');
end
assert(isnumeric(data)&&isreal(data)&&~issparse(data)&&ndims(data)<=3&& ...
    isequal([size(data,1) size(data,2) size(data,3)],double(op.projection_shape(:))'), ...
    'TIGRE:operator:shape','Invalid projection shape/type');
b=single(data); if ~isempty(mask), b(~mask)=0; end
assert(all(isfinite(b(:))),'TIGRE:operator:finite','Nonfinite valid measurement');
end
