% This file will compile all the necesary mex files for TIGRE to work. You
% need CUDA nvcc installed in your PC and setup with MATLAB mex.
%
%
% IMPORTANT:
% Due to hardware limitations we
% have been unable to test in all possible OS an dMATLAB versions.
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

cudapath=getenv('CUDA_PATH');


if isempty(cudapath)
    error(sprintf('CUDA Path not found. \nAdd the path by writting in MATLAB:\nsetenv(''CUDA_PATH'',''your path'')\nWhere "your path" is C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v8.0, for example')) ;
end

% Replace path fo current CUDA version in xml
% if ispc
%     fid  = fopen('mex_CUDA_win64.xml','r');
%     f=fread(fid,'*char')';
%     fclose(fid);
%     cudaverdiff=strfind(f,cudapath(end-2:end));
%     if ~isempty(cudaverdiff)
%         fid  = fopen('mex_CUDA_win64.xml','w');
%         f=strrep(f,'8.0',cudapath(end-2:end)); %the mex file has 8.0 on it
%         fwrite(fid,f);
%         fclose(fid);
%     end
%     
% end

% Compile for x64 or x32
disp('Compiling TIGRE source...')
disp('This may take a couple of minutes....')

if ispc
    % make sure the correc VS commons is in the file
%     msv10=getenv('VS100COMNTOOLS');
%     msv12=getenv('VS120COMNTOOLS');
%     msv14=getenv('VS140COMNTOOLS');
%     
%     if ~isempty(msv10)
%         fid  = fopen('mex_CUDA_win64.xml','r');
%         f=fread(fid,'*char')';
%         fclose(fid);
%         f=strrep(f,'VS120COMNTOOLS','VS100COMNTOOLS');
% 
%         fid  = fopen('mex_CUDA_win64.xml','w');
%         fwrite(fid,f);
%         fclose(fid);
%     end
%     if ~isempty(msv14)
%         fid  = fopen('mex_CUDA_win64.xml','r');
%         f=fread(fid,'*char')';
%         fclose(fid);
%         f=strrep(f,'VS120COMNTOOLS','VS140COMNTOOLS');
% 
%         fid  = fopen('mex_CUDA_win64.xml','w');
%         fwrite(fid,f);
%         fclose(fid);
%         
%         warning('You are using VS2015.')
%         disp('If you are using Visual Studio 2015 you may get errors in compilation.')
%         disp('If you the following error:')
%         disp('LINK : fatal error LNK1104: cannot open file ''ucrt.lib''')
%         disp('Do the following:')
%         disp(char(10))
%         disp('1: Locate your file ''ucrt.lib'' file, it would be somewhere similar to:')
%         disp('C:\Program Files (x86)\Windows Kits\10\Lib\10.0.10240.0\ucrt\x64      (The numbers may differ)')
%         disp(char(10))
%         disp('2: Open Compile.m')
%         disp(char(10))
%         disp('3: On lines 102-105, add in the end of the lines the following string with the correct path in your PC')
%         disp(' -I''C:\Program Files (x86)\Windows Kits\10\Include\10.0.10150.0\ucrt'' -L''C:\Program Files (x86)\Windows Kits\10\Lib\10.0.10150.0\ucrt\x64''')
%         disp(char(10))
%         disp('4: Try compiling again');
%         disp(char(10))
%         disp('If this does not work, copy ''ucrt.lib'' from the previous path to TIGRE folder and compile again (without the lines ons tep 3)')
%         disp(char(10))
%         disp('For any other error please contact the authors for help')
%         warning('END VS2015 WARNING')
%     end
%     
%     if(isempty(msv14)&&isempty(msv10)&&isempty(msv12))
%         error('VSCOMNTOOLS not found');
%     end
    
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\x64\cudart.lib" ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/win64
        mex -largeArrayDims "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\x64\cudart.lib" ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_spherical.cu ./Source/voxel_backprojection2_spherical.cu ./Source/voxel_backprojection_parallel.cu ./Source/voxel_backprojection_parallel_spherical.cu -outdir ./Mex_files/win64
        mex -largeArrayDims "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\x64\cudart.lib" ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/win64
        mex -largeArrayDims "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\x64\cudart.lib" ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/win64
        mex -largeArrayDims "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\x64\cudart.lib" ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/win64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/win64
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_spherical.cu ./Source/voxel_backprojection2_spherical.cu ./Source/voxel_backprojection_parallel.cu ./Source/voxel_backprojection_parallel_spherical.cu -outdir ./Mex_files/win64
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/win32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/win32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/win32
    end
    
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        disp('compiling for mac 64')
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_spherical.cu ./Source/voxel_backprojection2_spherical.cu ./Source/voxel_backprojection_parallel.cu ./Source/voxel_backprojection_parallel_spherical.cu -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/mac64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/mac64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/mac32
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_spherical.cu ./Source/voxel_backprojection2_spherical.cu ./Source/voxel_backprojection_parallel.cu ./Source/voxel_backprojection_parallel_spherical.cu -outdir ./Mex_files/mac32
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/mac32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/mac32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/mac32
    end
    
elseif isunix
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_spherical.cu ./Source/voxel_backprojection2_spherical.cu ./Source/voxel_backprojection_parallel.cu ./Source/voxel_backprojection_parallel_spherical.cu -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/linux64
        mex -largeArrayDims ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/linux64
    else
        mex  ./Source/Ax_mex.cpp ./Source/ray_interpolated_projection.cu ./Source/Siddon_projection.cu ./Source/ray_interpolated_projection_parallel.cu ./Source/Siddon_projection_parallel.cu -outdir ./Mex_files/linux32
        mex  ./Source/Atb_mex.cpp ./Source/voxel_backprojection.cu ./Source/voxel_backprojection2.cu ./Source/voxel_backprojection_spherical.cu ./Source/voxel_backprojection2_spherical.cu ./Source/voxel_backprojection_parallel.cu ./Source/voxel_backprojection_parallel_spherical.cu -outdir ./Mex_files/linux32
        mex  ./Source/minTV.cpp ./Source/POCS_TV.cu  -outdir ./Mex_files/linux32
        mex  ./Source/AwminTV.cpp ./Source/POCS_TV2.cu  -outdir ./Mex_files/linux32
        mex  ./Source/tvDenoise.cpp ./Source/tvdenoising.cu  -outdir ./Mex_files/linux32
    end
end





disp('')
disp('Compilation complete')
