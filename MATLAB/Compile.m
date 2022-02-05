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

if ispc 
    mex -setup:'mex_CUDA_win64.xml'
elseif ismac
    mex -setup:'mex_CUDA_maci64.xml'
elseif isunix    
    mex -setup:'mex_CUDA_glnxa64.xml'
end

addpath('./Utilities/Setup');

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

rmpath('./Utilities/Setup');

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
        mex -largeArrayDims ./Utilities/cuda_interface/AddNoise.cpp ../Common/CUDA/RandomNumberGenerator.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win64
    else
        mex  ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/win64
        mex  ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu  ../Common/CUDA/voxel_backprojection_parallel.cu -outdir ../Common/CUDA/GpuIds.cpp ./Mex_files/win64
        mex  ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/cuda_interface/AddNoise.cpp ../Common/CUDA/RandomNumberGenerator.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/IO/VarianCBCT/mexReadXim.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/win32
    end
    
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        disp('compiling for mac 64')
        mex -largeArrayDims ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu  ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/cuda_interface/AddNoise.cpp ../Common/CUDA/RandomNumberGenerator.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac64
    else
        mex  ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/cuda_interface/AddNoise.cpp ../Common/CUDA/RandomNumberGenerator.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/mac32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/mac32
    end
    
elseif isunix
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu  ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/cuda_interface/AddNoise.cpp ../Common/CUDA/RandomNumberGenerator.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux64
    else
        mex  ./Utilities/cuda_interface/Ax_mex.cpp ../Common/CUDA/ray_interpolated_projection.cu ../Common/CUDA/Siddon_projection.cu ../Common/CUDA/ray_interpolated_projection_parallel.cu ../Common/CUDA/Siddon_projection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/Atb_mex.cpp ../Common/CUDA/voxel_backprojection.cu ../Common/CUDA/voxel_backprojection2.cu ../Common/CUDA/voxel_backprojection_parallel.cu ../Common/CUDA/GpuIds.cpp -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/minTV.cpp ../Common/CUDA/POCS_TV.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/AwminTV.cpp ../Common/CUDA/POCS_TV2.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/tvDenoise.cpp ../Common/CUDA/tvdenoising.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/cuda_interface/AddNoise.cpp ../Common/CUDA/RandomNumberGenerator.cu ../Common/CUDA/GpuIds.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  -largeArrayDims ./Utilities/IO/VarianCBCT/mexReadXim.cpp -outdir ./Mex_files/linux32
        mex  ./Utilities/GPU/getGpuName_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
        mex  ./Utilities/GPU/getGpuCount_mex.cpp ../Common/CUDA/gpuUtils.cu -outdir ./Mex_files/linux32
    end
end

disp('')
disp('Compilation complete')

