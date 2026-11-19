function result = runIsolated(operation, backend, detail)
%RUNISOLATED Run a backend check in a separate MATLAB process.
%   Failed runs retain their request, result, and log in the reported folder.

if nargin < 3
    detail = '';
end
backend = validatestring(backend, {'cpu', 'metal', 'cuda'});
operation = validatestring(operation, ...
    {'clear', 'demo', 'recovery', 'fail', 'noResult'});
folder = tempname;
mkdir(folder);
requestFile = fullfile(folder, 'request.mat');
resultFile = fullfile(folder, 'result.mat');
logFile = fullfile(folder, 'matlab.log');
fprintf('Isolated %s/%s %s; log: %s\n', backend, operation, detail, logFile);
request = struct('operation', operation, 'backend', backend, ...
    'detail', detail, 'resultFile', resultFile, 'folder', folder);
save(requestFile, 'request');

testFolder = fileparts(fileparts(mfilename('fullpath')));
expression = sprintf("addpath('%s'); backendtest.worker('%s')", ...
    strrep(testFolder, "'", "''"), strrep(requestFile, "'", "''"));
executable = fullfile(matlabroot, 'bin', 'matlab');
if ispc
    executable = [executable, '.exe'];
end
argumentsList = java.util.ArrayList;
argumentsList.add(executable);
argumentsList.add('-batch');
argumentsList.add(char(expression));
builder = java.lang.ProcessBuilder(argumentsList);
builder.redirectErrorStream(true);
builder.redirectOutput(java.io.File(logFile));
process = builder.start();
processCleanup = onCleanup(@() stopProcess(process));
while process.isAlive()
    pause(0.25);
end
status = process.exitValue();
clear processCleanup

if status ~= 0
    error('TIGRE:Tests:ChildFailed', ...
        'Isolated %s check failed (exit %d). Log: %s\n%s', ...
        operation, status, logFile, fileread(logFile));
end
if ~isfile(resultFile)
    error('TIGRE:Tests:MissingChildResult', ...
        'Isolated %s check returned no result. Log: %s', operation, logFile);
end
loaded = load(resultFile, 'result');
if ~isfield(loaded, 'result') || ~loaded.result.passed
    error('TIGRE:Tests:InvalidChildResult', ...
        'Isolated %s check returned an invalid result. Log: %s', operation, logFile);
end
result = loaded.result;
rmdir(folder, 's');
end

function stopProcess(process)
if process.isAlive()
    process.destroyForcibly();
end
end
