function set_cuda_cc_flags(cuda_version)

if ispc
    filename='mex_CUDA_win64.xml';
elseif ismac
    filename='mex_CUDA_maci64.xml';
elseif isunix
    filename='mex_CUDA_glnxa64.xml';
end

fid=fopen(filename,'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

if ispc
    idx=squeeze(cellfun(@(x)strfind(x,'COMPILER="nvcc"'),A,'UniformOutput',false));
else % same for mac and linux
    idx=squeeze(cellfun(@(x)strfind(x,'MATLABMEX="-DMATLAB_MEX_FILE"'),A,'UniformOutput',false));
end
idx=find(~cellfun('isempty', idx));

if cuda_version>=11.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_86,code=sm_86 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>10.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 -gencode=arch=compute_75,code=sm_75 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>9.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>8.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>7.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
else
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
end
% Write cell A into txt
fid = fopen(filename, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end

end

