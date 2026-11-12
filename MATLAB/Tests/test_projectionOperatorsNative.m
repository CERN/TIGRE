function test_projectionOperatorsNative(outputFile)
% Native masked multi-panel rays, nested subsets, solver/wrapper parity.
if nargin<1, outputFile=''; end
geo=struct('nVoxel',[6;5;4],'dVoxel',[1.2;.8;1.5],'offOrigin',[.13;-.21;.17]);
beta=[.1 .7 1.6]; surface.source_positions_mm=zeros(3,3);
surface.pixel_centers_mm=zeros(3,5,3,3); surface.valid_mask=true(3,5,3);
[u,v]=meshgrid(linspace(-7,7,5),linspace(-4,4,3));
local=cat(3,-20+.03*u.^2,u,v);
for k=1:3
    R=[cos(beta(k)) -sin(beta(k)) 0;sin(beta(k)) cos(beta(k)) 0;0 0 1];
    surface.source_positions_mm(:,k)=R*[30;0;.3];
    surface.pixel_centers_mm(:,:,:,k)=reshape(reshape(local,[],3)*R',3,5,3);
end
surface.valid_mask(:,3,:)=false;
op=projectionOperator(geo,surface,'Device',0);
x=zeros(geo.nVoxel','single'); x(2:5,2:4,2:3)=1; x(4,3,3)=1.7;
data=op.forward(x); data(~op.valid_mask)=NaN;
sub=op.subset([3 1 3]); nested=sub.subset([2 1]);
full=op.forward(x); selected=nested.forward(x); expected=full(:,:,[1 3]);
assert(max(abs(selected(:)-expected(:)))<2e-5);
q=reshape(single(sin(1:numel(data))),size(data)); q(~op.valid_mask)=0;
at=op.adjoint(q); left=sum(double(full(:)).*double(q(:))); right=sum(double(x(:)).*double(at(:)));
assert(abs(left-right)<2e-5*max([1 abs(left) abs(right)]));
names={'CGLS', 'SIRT', 'SART', 'OS_SART', 'fast_OS_SART', 'SART_TV', 'OS_SART_TV', 'iterativereconalg'}; results=struct;
for item=names
    name=item{1}; options={};
    if contains(lower(name),'pocs')||contains(lower(name),'pcsd')
        options={'MaxL2err',0,'TViter',2};
    elseif contains(lower(name),'sart_tv')
        options={'TViter',2,'TVLambda',50};
    elseif ismember(name,{'FISTA','ISTA'})
        options={'Hyper',200,'TViter',2,'lambda',.1};
    elseif strcmp(name,'IRN_TV_CGLS')
        options={'niter_outer',1};
    end
    actual=feval([name '_operator'],data,op,2,options{:});
    if exist(name,'file')==2
        public=feval(name,data,[],[],2,'Operator',op,options{:});
        assert(norm(double(actual(:))-double(public(:)))<2e-3*max(norm(double(actual(:))),1),'Public mismatch %s',name);
    end
    wrapped=feval([name '_surface'],data,geo,surface,2,options{:},'Device',0);
    assert(isa(actual,'single')&&all(isfinite(actual(:))),'Nonfinite %s',name);
    assert(norm(double(actual(:))-double(wrapped(:)))<2e-3*max(norm(double(actual(:))),1),'Wrapper mismatch %s',name);
    residual=op.forward(actual)-operatorPrepareData(data,op);
    assert(norm(double(residual(:)))<1.1*norm(double(full(:))),'No progress %s',name);
    results.(name)=actual;
end
if ~isempty(outputFile), save(outputFile,'geo','surface','data','results','-v7'); end
fprintf('PASS native operator/wrapper parity: %d methods; MEX %s\n',numel(names),which('Surface_mex'));
end
