function test_projectionOperators()
% Independent matrix references with a missing ray and unobservable voxel.
A=single([1 2 0 0;0 1 3 0;2 0 1 0;1 1 1 0;0 2 1 0;3 1 0 0]); A(3,:)=0;
op=matrixOp(A); truth=reshape(single([1 2 3 0]),[2 2 1]); b=op.forward(truth);
b(~op.valid_mask)=NaN;
[x,history]=CGLS_operator(b,op,3);
expected=pinv(double(A))*reshape(double(operatorPrepareData(b,op)),[],1);
assert(norm(double(x(:))-expected(:))<2e-5); assert(numel(history)<=3);
clean=operatorPrepareData(b,op); row=sum(double(A),2); col=sum(double(A),1)';
r=zeros(6,1); r(row>0)=double(clean(row>0))./row(row>0);
expected=zeros(4,1); correction=double(A)'*r; expected(col>0)=correction(col>0)./col(col>0);
x=SIRT_operator(b,op,1); assert(norm(double(x(:))-expected)<2e-6);
init=ones(2,2,'single'); saved=init;
x=CGLS_operator(op.forward(init),op,3,'Init',init);
assert(isequal(x,init)&&isequal(saved,init));
for bad={-1,1.5}
    mustFail(@() CGLS_operator(b,op,bad{1}));
end
mustFail(@() CGLS_operator(b,op,0,'UnknownOption',1));
mustFail(@() CGLS_operator(b,op,0,'Init',nan(2,2)));
sub=op.subset([3 1 3]); nested=sub.subset([2 1]);
q=op.forward(truth); z=nested.forward(truth); assert(isequal(z,q(:,:,[1 3])));
names={'CGLS', 'SIRT', 'SART', 'OS_SART', 'fast_OS_SART', 'SART_TV', 'OS_SART_TV', 'iterativereconalg'};
for item=names
    name=item{1}; options={}; count=2;
    if contains(lower(name),'pocs')||contains(lower(name),'pcsd')
        options={'MaxL2err',0,'TViter',0};
    elseif contains(lower(name),'sart_tv')
        options={'TViter',0};
    elseif ismember(name,{'FISTA','ISTA'})
        options={'Hyper',200,'TViter',0,'lambda',0};
    elseif strcmp(name,'IRN_TV_CGLS')
        count=0; % Native regularizer requires all volume axes >= 2.
    elseif strcmp(name,'hybrid_fLSQR_TV')
        count=0;
    end
    x=feval([name '_operator'],b,op,count,options{:});
    assert(isa(x,'single')&&all(isfinite(x(:))),'Nonfinite %s',name);
    if count>0
        r=op.forward(x)-clean; assert(norm(double(r(:)))<1.1*norm(double(clean(:))),'No progress %s',name);
    end
end
fprintf('PASS shared operators: dense oracles, masks, subsets, validation; %d methods\n',numel(names));
end
function op=matrixOp(A)
op.volume_shape=[2 2 1]; op.projection_shape=[1 2 size(A,1)/2];
op.valid_mask=reshape(any(A,2),op.projection_shape);
op.forward=@(x) reshape(A*x(:),op.projection_shape);
op.adjoint=@(b) reshape(A'*reshape(single(b).*single(op.valid_mask),[],1),op.volume_shape);
op.subset=@(idx) matrixSubset(A,idx);
end
function op=matrixSubset(A,idx)
rows=reshape(1:size(A,1),2,[]); rows=rows(:,idx); op=matrixOp(A(rows(:),:));
end
function mustFail(action)
failed=false; try, action(); catch, failed=true; end
assert(failed,'Invalid input must fail');
end
