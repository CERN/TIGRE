function set_cuda_cc_flags(cuda_version, cuda_path)
% Write the nvcc flags, including the GPU architectures to compile for, into
% the mex_CUDA_*.xml file that mex uses.
%
% The architectures are the ones this nvcc supports, as it reports them
% (nvcc --list-gpu-arch), the same as the Python build does (setup.py). The
% fixed table further down is only for an nvcc too old to report them, or
% when cuda_path is not given.

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

gencode='';
if nargin>1 && ~isempty(cuda_path)
    gencode=gencode_from_nvcc(cuda_path);
end
if ~isempty(gencode)
   A{idx+1}=['        COMPFLAGS=" ', gencode, ' --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"'];
elseif cuda_version>=13.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_89,code=sm_89 -gencode=arch=compute_90,code=sm_90 -gencode=arch=compute_120,code=sm_120 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>=12.8  % compute_120 needs CUDA 12.8
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_89,code=sm_89 -gencode=arch=compute_90,code=sm_90 -gencode=arch=compute_120,code=sm_120 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';  
elseif cuda_version>=12.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_89,code=sm_89 -gencode=arch=compute_90,code=sm_90 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>=11.8
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_89,code=sm_89 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>=11.1
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
elseif cuda_version>=11.0
   A{idx+1}='        COMPFLAGS=" -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"';
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

function gencode=gencode_from_nvcc(cuda_path)
% -gencode flags for every architecture this nvcc supports, or '' if it
% cannot tell (nvcc older than --list-gpu-arch).
gencode='';
if ispc
    [status,cmout]=system(['"', cuda_path, '/bin/nvcc" --list-gpu-arch']);
else
    [status,cmout]=system([cuda_path, '/bin/nvcc --list-gpu-arch']);
end
if status
    return
end
% Plain architectures only (compute_90, not the specific compute_90a/100f).
tok=regexp(cmout,'^compute_(\d+)\s*$','tokens','lineanchors');
if isempty(tok)
    return
end
cc=unique(cellfun(@(t) str2double(t{1}), tok));
for c=cc
    gencode=[gencode, sprintf('-gencode=arch=compute_%d,code=sm_%d ', c, c)]; %#ok<AGROW>
end
gencode=strtrim(gencode);
end
