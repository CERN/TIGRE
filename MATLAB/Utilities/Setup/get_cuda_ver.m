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
