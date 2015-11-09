function [ data ] = ParkerWeight( data, param )
%PARKERWEIGHT Summary of this function goes here
%   Detailed explanation goes here

as = atan(param.us/param.DSD)/pi*180;

fanang = abs(as(end)-as(1))/2;

for  i=1:param.nProj
    pweight  = ((180/param.dang - param.nProj)/2+i)/(180/param.dang)*180 + param.dir*as;
    
    id  = find(pweight<=fanang/2 & pweight>-fanang/2); 
    if  isempty(id)  ~= 1
        data(id,:,i) = data(id,:,i).*repmat(max(sin(pi/2*(fanang/2+pweight(id)')./fanang),0).^2,[1, param.nv,1]);
    end
    
    id = find(pweight>=180-fanang/2 & pweight<180+fanang/2);
    if  isempty(id)  ~= 1
         data(id,:,i) = data(id,:,i).*repmat(max(cos(pi/2*(fanang/2+pweight(id)'-180)./fanang),0).^2,[1, param.nv,1]);
    end
    
    id  = find(pweight<-fanang/2 | pweight>180+fanang/2); % Exception check
    if  isempty(id)  ~= 1
        data(id,:,i) = 0;
    end
    
end

end

