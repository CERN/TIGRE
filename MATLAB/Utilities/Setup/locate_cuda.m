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
