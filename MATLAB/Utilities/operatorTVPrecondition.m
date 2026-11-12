function x=operatorTVPrecondition(weight,rhs)
% Bounded CG inversion of D' W^2 D on its mean-free subspace.
rhs=single(reshape(rhs,size(weight))); rhs=rhs-single(mean(double(rhs(:))));
x=zeros(size(rhs),'single'); r=rhs; p=r;
gamma=squared(r); initial=gamma;
for iteration=1:50
    if gamma==0 || gamma<=1e-12*initial, break; end
    q=zeros(size(rhs),'single');
    for axis=1:3
        order=[axis setdiff(1:3,axis)]; value=permute(p,order); w=permute(weight,order);
        gradient=zeros(size(value),'single'); gradient(1:end-1,:,:)=value(1:end-1,:,:)-value(2:end,:,:);
        value=w.^2.*gradient; adjoint=zeros(size(value),'single');
        adjoint(1:end-1,:,:)=adjoint(1:end-1,:,:)+value(1:end-1,:,:);
        adjoint(2:end,:,:)=adjoint(2:end,:,:)-value(1:end-1,:,:);
        q=q+ipermute(adjoint,order);
    end
    denominator=sum(double(p(:)).*double(q(:))); if denominator<=0, break; end
    alpha=gamma/denominator; x=x+single(alpha)*p; r=r-single(alpha)*q;
    next=squared(r); p=r+single(next/gamma)*p; gamma=next;
end
end
function result=squared(x)
result=sum(double(x(:)).^2);
end
