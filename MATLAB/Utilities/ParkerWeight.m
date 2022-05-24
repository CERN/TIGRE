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
if length(unique(geo.DSD))>1
   warning('Parker weigths not supported for varying geo.DSD, mean(geo.DSD) is used.') % mean(DSD) might be better than DSD(1) for varying DSD
end

alpha = atan([-geo.sDetector(1)/2+geo.dDetector(1)/2:geo.dDetector(1):geo.sDetector(1)/2-geo.dDetector(1)/2]/mean(geo.DSD));
diff_angles = diff(angles);
if (all(diff_angles>=0)) % if angles are incremental, flip alpha in the left-right direction
   alpha=-alpha;
end
delta = abs(alpha(end)-alpha(1))/2;
totangles=abs(angles(end)-angles(1));
if totangles>=2*pi
   warning('Computing Parker weigths for scanning angle equal or bigger than 2*pi. Consider disabling Parker weigths.') 
end
if totangles<pi+2*delta
    warning('Scanning angles smaller than pi+cone_angle. This is limited angle tomgraphy, there is no sufficient data, thus weighting for data redundancy is not required.')
end
epsilon=max(totangles-(pi+2*delta),0);



for ii=1:size(data,3)
    beta=abs(angles(ii)-angles(1)); % make sure beta always starts from 0 to totangles
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
