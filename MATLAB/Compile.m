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
    cudapath=getenv('CUDA_PATH');
    if isempty(cudapath)
        error(sprintf('CUDA Path not found. \nAdd the path by writting in MATLAB:\nsetenv(''CUDA_PATH'',''your path'')\nWhere "your path" is C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v8.0, for example')) ;
    end
end

% Compile for x64 or x32
disp('Compiling TIGRE source...')
disp('This may take a couple of minutes....')

if ispc

    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/win64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/win64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/win64
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu -outdir ./Mex_files/win64
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/win32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/win32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/win32
    end
    
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        disp('compiling for mac 64')
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu  -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/mac64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/mac32
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_parallel.cu  -outdir ./Mex_files/mac32
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/mac32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/mac32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/mac32
    end
    
elseif isunix
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu  ./Source/voxel_backprojection_parallel.cu  -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/linux64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/linux32
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_parallel.cu  -outdir ./Mex_files/linux32
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/linux32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/linux32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/linux32
    end
end





disp('')
disp('Compilation complete')
