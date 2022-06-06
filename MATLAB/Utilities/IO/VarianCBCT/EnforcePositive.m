function proj = EnforcePositive(proj)
%% Remove anomalies
% in case of NaN or Inf
proj(isnan(proj)) = 0;
proj(isinf(proj)) = 0;

%% Set all negative to zeros
proj(proj<0) = 0;

end
