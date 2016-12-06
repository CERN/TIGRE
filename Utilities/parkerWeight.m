function [ data ] = parkerWeight( data, geo,angles )
%PARKERWEIGHT Summary of this function goes here
%   Detailed explanation goes here

as = atan([-geo.sDetector(1)/2+geo.dDetector(1)/2:geo.dDetector(1):geo.sDetector(1)/2-geo.dDetector(1)/2]/geo.DSD)/pi*180;
fanang = abs(as(end)-as(1))/2;

dang=diff(angles*180/pi);
ndangle=(unique(dang));
if any(diff(ndangle)>1e-10)
    warning('Parker wegths assume uniform angles, and these ones dont seem that they are');
end
dang=dang(1);

for  ii=1:size(data,3)
    pweight  = ((180/abs(dang) - size(data,3))/2+ii)/(180/abs(dang))*180 - sign(dang)*as;
    
    id  = find(pweight<=fanang/2 & pweight>-fanang/2); 
    if  isempty(id)  ~= 1
        data(:,id,ii) = data(:,id,ii).*repmat(max(sin(pi/2*(fanang/2+pweight(id)')./fanang),0).^2,[1, geo.nDetector(2),1]).';
    end
    
    id = find(pweight>=180-fanang/2 & pweight<180+fanang/2);
    if  isempty(id)  ~= 1
         data(:,id,ii) = data(:,id,ii).*repmat(max(cos(pi/2*(fanang/2+pweight(id)'-180)./fanang),0).^2,[1, geo.nDetector(2),1]).';
    end
    
    id  = find(pweight<-fanang/2 | pweight>180+fanang/2); % Exception check
    if  isempty(id)  ~= 1
        data(:,id,ii) = 0;
    end
    
end

end

