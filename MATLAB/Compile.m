% This file will compile all the necesary mex files for TIGRE to work. You
% need CUDA nvcc installed in your PC and setup with MATLAB mex.
%
%
% IMPORTANT:
% Due to hardware limitations we
% have been unable to test in all possible OS and MATLAB versions.
% Please, report any issue with compilation in other systems
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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
%% Clear all clears also mex
clear all;
%% Compile

mex -setup

if ispc
    currentFolder = cd;
    fileExisting  = (exist(fullfile(currentFolder, 'mex_CUDA_win64.xml'), 'file') == 2);
    if ~fileExisting
        error(sprintf('mex_CUDA_win64.xml not found. You may need to rename the existing files depending on your MVS version')) ;
    end
end
[cudapath, cuda_ver]=locate_cuda();
if isempty(cudapath)
    error(sprintf('CUDA Path not found. \nAdd the path by writting in MATLAB:\nsetenv(''CUDA_PATH'',''your path'')\nWhere "your path" is C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.2, for example, \nor /usr/local/cuda on linux')) ;
end
if ispc
   setenv('CUDA_PATH',cudapath);
end

set_cuda_cc_flags(cuda_ver);

% Compile for x64 or x32
disp('Compiling TIGRE source...')
disp('This may take a couple of minutes....')

if ispc
    
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu  ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/cuda_interface/minPICCS.cpp ../Common/CUDA/PICCS.cu  -outdir ./Mex_files/win64

    else
        mex  ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/win64
        mex  ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu  ../Common/CUDA/voxel_backprojection_parallel.cu -outdir ../Common/CUDA/GpuIds.cpp ./Mex_files/win64
        mex  ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/win32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/cuda_interface/minPICCS.cpp ../Common/CUDA/PICCS.cu  -outdir ./Mex_files/win32

    end
    
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        disp('compiling for mac 64')
        mex -largeArrayDims ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu  ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/minPICCS.cpp ../Common/CUDA/PICCS.cu  -outdir ./Mex_files/mac64
    else
        mex  ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/mac32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/minPICCS.cpp ../Common/CUDA/PICCS.cu  -outdir ./Mex_files/win32

    end
    
elseif isunix
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu  ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/minPICCS.cpp ../Common/CUDA/PICCS.cu  -outdir ./Mex_files/linux64

    else
        mex  ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/linux32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/minPICCS.cpp ../Common/CUDA/PICCS.cu  -outdir ./Mex_files/linux32
    end
end

disp('')
disp('Compilation complete')


function [cuda_path, cuda_ver]=locate_cuda()

cuda_ver=-1;
% Guess 1:
cuda_path=getenv('CUDA_PATH');
if isempty(cuda_path)
    cuda_path=getenv('CUDA_HOME');
end
if ~isempty(cuda_path) % we have something.
    cuda_ver=get_cuda_ver(cuda_path);
    return
end
% Guess 2:
if ispc
    which='where';
else
    which='which';
end
[status,cmout]=system([which, ' nvcc']);
if ~status % succeded
    verstr=strsplit(cmout,'\n');
    %which one to use? the first one I guess.
    verstr=verstr{1};
    cuda_path=strsplit(verstr,'bin');  
    cuda_path=cuda_path{1}(1:end-1);
    cuda_ver=get_cuda_ver(cuda_path);
    return
end
% Guess 3
if ispc
    guess_cuda_path='C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/';
    if exist(guess_cuda_path, 'dir')
        versions=ls(guess_cuda_path);
        % just grab biggest one.
        versions_num=str2double(versions(3:end,2:end));
        [~,idx]=max(versions_num);
        cuda_path=[guess_cuda_path,versions(2+idx,:)];
        cuda_ver=get_cuda_ver(cuda_path);
        return
    end
else
    % symlinc
    guess_cuda_path='/usr/local/cuda';
    if exist(guess_cuda_path, 'dir')
        cuda_path=guess_cuda_path;
        cuda_ver=get_cuda_ver(cuda_path);
        return
    end
end


end
function cuda_ver=get_cuda_ver(cuda_path)
if ispc
    [status,cmout]=system(['"', cuda_path, '/bin/nvcc" -V']);
else
    [status,cmout]=system([cuda_path, '/bin/nvcc -V']);
end
if status
    error('Error finding CUDA version')
else
    stridx=strfind(cmout,'release ');
    cuda_ver=str2double(cmout(stridx+length('release ') : stridx+length('release ')+3));
end
end
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