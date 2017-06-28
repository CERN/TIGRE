function  citeme( func_name ,bib)
%CITEME(FUCTIONAME) gives the bilbiografic reference of the fucntion being
% used.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
if nargin<2
    bib=false;
end
disp('')

file=fopen('TIGRE.bib');
tline = fgets(file);
found=false;
while ischar(tline)&&~found
    tline=strtrim(tline);
    if strcmpi(tline,func_name)
        found=true;
        break;
    end
    tline = fgetl(file);
end
% found the bib entry

if found
    tline = fgets(file);
    if ~bib
        while  ~strcmpi(tline(1),'@')
            disp(tline(1:end-1));
            tline = fgets(file);
        end
    else
        while ~strcmpi(tline(1),'@')
            tline = fgets(file);
        end
        while ~isempty(strtrim(tline))
            disp(tline(1:end-1));
            tline = fgets(file);
        end
    end
    
else
   disp('The function needs no bibliographic reference (other than TIGRE)')
   disp([' [' 8 ' Warning: Make sure there is no spelling mistake.]' 8 ''])
end
fclose(file);
disp('')

end

