function backend = tigreBackend()
%TIGREBACKEND Return the active TIGRE compute backend.
%   BACKEND = tigreBackend() returns "cuda", "cpu", or "metal".

if ~isappdata(groot, 'TIGREBackendState')
    error('TIGRE:Backend:NotInitialized', ...
        'TIGRE is not initialized. Run InitTIGRE before using the toolbox.');
end

state = getappdata(groot, 'TIGREBackendState');
backend = string(state.backend);
end
