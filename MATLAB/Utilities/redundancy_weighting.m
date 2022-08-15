function w = redundancy_weighting(geo)
% Preweighting using Wang function
% Ref: 
%    Wang, Ge. X-ray micro-CT with a displaced detector array. Medical Physics, 2002,29(7):1634-1636.

if ~isfield(geo,'COR')
    geo.COR=0;
end
offset = geo.offDetector(1);
offset = offset + (geo.DSD(1) / geo.DSO(1)) * geo.COR(1);   % added correction
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + abs(offset);

us = us * geo.DSO(1)/geo.DSD(1);
theta = (geo.sDetector(1)/2 - abs(offset))...
        * sign(offset);
abstheta = abs(theta * geo.DSO(1)/geo.DSD(1));

w = ones([geo.nDetector(2),geo.nDetector(1)]);
if apply_wang_weights(geo)
    for ii = 1:geo.nDetector
        t = us(ii);
        if(abs(t) <= abstheta)
            w(:,ii) = 0.5*(sin((pi/2)*atan(t/geo.DSO(1))/(atan(abstheta/geo.DSO(1)))) + 1);
        end
        if(t<-abstheta)
            w(:,ii) = 0;
        end
    end
    w=w*2;
    if(theta<0)
        w = fliplr(w);
    end
end
end