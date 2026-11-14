function output=operatorProject(back,data,op,angles)
% Internal view tokens are indices, never physical angles for native kernels.
idx=double(angles(1,:));
validateattributes(idx,{'double'},{'integer','>=',1,'<=',op.projection_shape(3),'finite'});
if isequal(idx,1:op.projection_shape(3)), sub=op; else, sub=op.subset(idx); end
if back
    data=single(data); if ~isempty(sub.valid_mask), data(~sub.valid_mask)=0; end
    assert(all(isfinite(data(:))),'TIGRE:operator:finite','Nonfinite valid data');
    output=single(sub.adjoint(data));
else
    output=single(sub.forward(single(data)));
    if ~isempty(sub.valid_mask), output(~sub.valid_mask)=0; end
end
assert(all(isfinite(output(:))),'TIGRE:operator:finite','Nonfinite operator result');
end
