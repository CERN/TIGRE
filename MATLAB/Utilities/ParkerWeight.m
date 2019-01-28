function data=ParkerWeight(data,geo,angles,q)


%Wesarg, Stefan, Matthias Ebert, and Thomas Bortfeld. 
% "Parker weights revisited." Medical physics 29.3 (2002): 372-378.

% @article{wesarg2002parker,
%   title={Parker weights revisited},
%   author={Wesarg, Stefan and Ebert, Matthias and Bortfeld, Thomas},
%   journal={Medical physics},
%   volume={29},
%   number={3},
%   pages={372--378},
%   year={2002},
%   publisher={American Association of Physicists in Medicine}
% }
if unique(geo.DSD)>1
   warning('Parker weigths not supported for varying geo.DSD, first one used.') 
end

alpha = atan([-geo.sDetector(1)/2+geo.dDetector(1)/2:geo.dDetector(1):geo.sDetector(1)/2-geo.dDetector(1)/2]/geo.DSD(1));
alpha=-alpha;
delta = abs(alpha(end)-alpha(1))/2;
totangles=cumsum(diff(angles));
totangles=totangles(end);
if totangles>=2*pi
   warning('Computing Parker weigths for scanning angle equal or bigger than 2*pi. Consider disabling Parker weigths.') 
end
if totangles<pi+2*delta
    warning('Scanning angles smaller than pi+cone_angle. This is limited angle tomgraphy, there is nosufficient data, thus weigthing for data redundancy is not required.')
end
epsilon=max(totangles-(pi+2*delta),0);



for ii=1:size(data,3)
    beta=angles(ii);
     w=0.5*(S(beta./b(alpha,delta,epsilon,q)-0.5)+S((beta-2*delta+2*alpha-epsilon)./b(alpha,delta,epsilon,q)+0.5)...
         -S((beta-pi+2*alpha)./b(-alpha,delta,epsilon,q)-0.5) ...
          -S((beta-pi-2*delta-epsilon)./b(-alpha,delta,epsilon,q)+0.5)...
          );
     data(:,:,ii)=data(:,:,ii).*repmat(w,[size(data,1),1]);
end

end

function w=S(beta)
w=zeros(size(beta));
w(beta<=-0.5)=0;
w(abs(beta)<0.5)=0.5*(1+sin(pi*beta(abs(beta)<0.5)));
w(beta>=0.5)=1;
end
function res=B(alpha,delta,epsilon)
    res=2*delta-2*alpha+epsilon;
end
function res=b(alpha,delta,epsilon,q)
    res=q*B(alpha,delta,epsilon);
end