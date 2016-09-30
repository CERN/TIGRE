function isit=is2014bOrNewer()
vers=version('-release');
if str2num(vers(1:4))>2014 || str2num(vers(1:4))==2014 && strcmp(vers(5),'b')
    isit=true;
else
    isit=false;
end

end