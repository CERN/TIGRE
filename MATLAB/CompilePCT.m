%Compile pCT radiography

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
% Coded by:           Stefanie Kaser, Benjamin Kirchmayer 
%--------------------------------------------------------------------------
clear all;

%% Precompile actions to write intercepts vector lengths
prompt = 'Enter pixel side length in mm (or hit enter): ';
input0 = input(prompt);


if input0 > 0
    pix = single(input0/0.25);

    prompt = 'Enter SL intercepts before phantom (or hit enter for default values): ';
    input1 = input(prompt);
    prompt = 'Enter CS intercepts within hull (or hit enter for default values): ';
    input2 = input(prompt);
    prompt = 'Enter SL intercepts behind phantom (or hit enter for default values): ';
    input3 = input(prompt);

    try 
        if input1 > 0
            interceps_in = uint16(input1);
        else
            interceps_in = uint16(10/pix);
            disp(strcat('Setting intercepts vector size (SL, in) to default: ', int2str(uint16(10/pix))));
        end
    catch
        interceps_in = uint16(10/pix);
        disp(strcat('Setting intercepts vector size (SL, in) to default: ', int2str(uint16(10/pix))));
    end

    try 
        if input2 > 0
            interceps_cs = uint16(input2);
        else
            interceps_cs = uint16(220/pix);
            disp(strcat('Setting intercepts vector size (CS within hull) to default: ', int2str(uint16(220/pix))));
        end
    catch
        interceps_cs = uint16(220/pix);
        disp(strcat('Setting intercepts vector size (CS within hull) to default: ', int2str(uint16(220/pix))));
    end

    try   
        if input3 > 0
            interceps_out = uint16(input3);
        else
            interceps_out = 100;
            disp(strcat('Setting intercepts vector size (SL, out) to default: ', int2str(uint16(100/pix))));
        end
    catch
        interceps_out = 100;
        disp(strcat('Setting intercepts vector size (SL, out) to default: ', int2str(uint16(100/pix))));
    end
else
    warning('No valid pixel size was entered. Intercepts vector size can be set manually or default values can be used. Be aware that default values are optimal for a pixel size of 0.25x0.25mm^2.')
    prompt = 'Enter SL intercepts before phantom (or hit enter for default values): ';
    input1 = input(prompt);
    prompt = 'Enter CS intercepts within hull (or hit enter for default values): ';
    input2 = input(prompt);
    prompt = 'Enter SL intercepts behind phantom (or hit enter for default values): ';
    input3 = input(prompt);

    try 
        if input1 > 0
            interceps_in = uint16(input1);
        else
            interceps_in = 10;
            disp('Setting intercepts vector size (SL, in) to default (10)');
        end
    catch
        interceps_in = 10;
        disp('Setting intercepts vector size (SL, in) to default (10)');
    end

    try 
        if input2 > 0
            interceps_cs = uint16(input2);
        else
            interceps_cs = 220;
            disp('Setting intercepts vector size (CS within hull) to default (220).');
        end
    catch
        interceps_cs = 220;
        disp('Setting intercepts vector size (CS within hull) to default (220).');
    end

    try   
        if input3 > 0
            interceps_out = uint16(input3);
        else
            interceps_out = 100;
            disp('Setting intercepts vector size (SL, out) to default (100).');
        end
    catch
        interceps_out = 100;
        disp('Setting intercepts vector size (SL, out) to default (100).');
    end
end

% SL (before phantom)
try
    textCell = readlines('../Common/CUDA/improvedForwardProjections.hpp');
catch
   error('Assure to start the compile script from the correct directory') 
end
searchMask = cell2mat(cellfun(@(x) contains(x, '#define vecSizeIn'), textCell','UniformOutput', false))';
textCell{searchMask} = ['#define vecSizeIn ' num2str(interceps_in)];

writeID = fopen('../Common/CUDA/improvedForwardProjections.hpp', 'w');

for i=1:numel(textCell)
    fprintf(writeID, '%s\n', textCell{i});
end
fclose(writeID);


% CS (within hull)
textCell = readlines('../Common/CUDA/improvedForwardProjections.hpp');
searchMask = cell2mat(cellfun(@(x) contains(x, '#define vecSizeCS'), textCell','UniformOutput', false))';
textCell{searchMask} = ['#define vecSizeCS ' num2str(interceps_cs)];

writeID = fopen('../Common/CUDA/improvedForwardProjections.hpp', 'w');

for i=1:numel(textCell)
    fprintf(writeID, '%s\n', textCell{i});
end
fclose(writeID);

% SL (behind phantom)
textCell = readlines('../Common/CUDA/improvedForwardProjections.hpp');
searchMask = cell2mat(cellfun(@(x) contains(x, '#define vecSizeOut'), textCell','UniformOutput', false))';
textCell{searchMask} = ['#define vecSizeOut ' num2str(interceps_out)];

writeID = fopen('../Common/CUDA/improvedForwardProjections.hpp', 'w');

for i=1:numel(textCell)
    fprintf(writeID, '%s\n', textCell{i});
end
fclose(writeID);


%% Compile
clear all;
addpath('pCTMexFiles');
addpath('./Utilities/Setup');
mex -setup
[cudapath, cuda_ver]=locate_cuda();
if isempty(cudapath)
    error(sprintf('CUDA Path not found. \nAdd the path by writting in MATLAB:\nsetenv(''CUDA_PATH'',''your path'')\nWhere "your path" is C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.2, for example, \nor /usr/local/cuda on linux')) ;
end
if ispc
   setenv('CUDA_PATH',cudapath);
end
set_cuda_cc_flags(cuda_ver);

rmpath('./Utilities/Setup');
% if isempty(getenv('CUDA_PATH'))
%    setenv('CUDA_PATH','usr/local/cuda') ;
% end
 
fprintf("Compiling pCT source...\n");
if isunix
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Utilities/cuda_interface/pCTCubicSpline_mex.cpp ../Common/CUDA/improvedForwardProjections.cu ../Common/CUDA/improvedForwardProjections_cone.cu -outdir ./pCTMexFiles/linux64
    else
        mex ./Utilities/cuda_interface/pCTCubicSpline_mex.cpp ../Common/CUDA/improvedForwardProjections.cu ../Common/CUDA/improvedForwardProjections_cone.cu -outdir ./pCTMexFiles/linux32
    end
elseif ispc
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Utilities/cuda_interface/pCTCubicSpline_mex.cpp ../Common/CUDA/improvedForwardProjections.cu ../Common/CUDA/improvedForwardProjections_cone.cu -outdir ./pCTMexFiles/win64
    else
        mex ./Utilities/cuda_interface/pCTCubicSpline_mex.cpp ../Common/CUDA/improvedForwardProjections.cu ../Common/CUDA/improvedForwardProjections_cone.cu -outdir ./pCTMexFiles/win32
    end
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        mex -largeArrayDims ./Utilities/cuda_interface/pCTCubicSpline_mex.cpp ../Common/CUDA/improvedForwardProjections.cu ../Common/CUDA/improvedForwardProjections_cone.cu -outdir ./pCTMexFiles/mac64
    else
        mex ./Utilities/cuda_interface/pCTCubicSpline_mex.cpp ../Common/CUDA/improvedForwardProjections.cu ../Common/CUDA/improvedForwardProjections_cone.cu -outdir ./pCTMexFiles/mac32
    end
end

if ispc
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./pCTMexFiles/win64');
    else
        addpath('./pCTMexFiles/win32');
    end
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./pCTMexFiles/mac64');
    else
        addpath('./pCTMexFiles/mac32');
    end
else
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./pCTMexFiles/linux64');
    else
        addpath('./pCTMexFiles/linux32');
    end
end

fprintf('Compilation successful!\n');