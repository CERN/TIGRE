function km=poissonrandom(lambda)

% Donald Knuth's poisson random generator

% This sucks bananas because its O(lambda) and we want lambda~=10000...

km=zeros(size(lambda));
for ii=1:numel(lambda)
    
    L=exp(-lambda(ii));
    k=1;
    p=1*rand;
    
    while p>L
        k=k+1;
        p=p*rand;
    end
    km(ii)=k-1;
end
end