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

% Compile for x64 or x32
disp('Compiling TIGRE source...')
disp('This may take a couple of minutes....')

if ispc
    
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/win64
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu -outdir ./Source/GpuIds.cpp ./Mex_files/win64
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/win32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/win32
    end
    
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        disp('compiling for mac 64')
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/mac32
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/mac32
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/mac32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/mac32
    end
    
elseif isunix
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/linux32
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_parallel.cu ./Source/GpuIds.cpp -outdir ./Mex_files/linux32
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu ./Source/GpuIds.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/linux32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ./Source/gpuUtils.cu -outdir ./Mex_files/linux32
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