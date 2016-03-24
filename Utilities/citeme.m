function  citeme( func_name ,bib)
%CITEME(FUCTIONAME) gives the bilbiografic reference of the fucntion being
% used.
if nargin<2
    bib=false;
end
disp('')

file=fopen('TIGRE.bib');
tline = fgets(file);
found=false;
while ischar(tline)&&~found
    if strcmpi(tline(1:end-1),func_name)
        found=true;
        break;
    end
    tline = fgets(file);
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
   disp('Warning: Make sure there is no spelling mistake.')
end
fclose(file);
disp('')

end

