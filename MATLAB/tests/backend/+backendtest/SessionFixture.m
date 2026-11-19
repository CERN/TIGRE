classdef SessionFixture < matlab.unittest.fixtures.Fixture
    %SessionFixture Initialise TIGRE and restore the caller's session.

    methods
        function setup(fixture)
            previousPath = path;
            previousFolder = pwd;
            previousRng = rng;
            previousVisibility = get(groot, 'DefaultFigureVisible');
            previousFigures = findall(groot, 'Type', 'figure');
            hadBackend = isappdata(groot, 'TIGREBackendState');
            previousBackend = [];
            if hadBackend
                previousBackend = getappdata(groot, 'TIGREBackendState');
            end
            fixture.addTeardown(@() backendtest.SessionFixture.restore( ...
                previousPath, previousFolder, previousRng, ...
                previousVisibility, previousFigures, hadBackend, previousBackend));

            matlabRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            addpath(matlabRoot, fullfile(matlabRoot, 'tests', 'fixtures'));
            InitTIGRE('backend', 'cpu');
            rng(42, 'twister');
            set(groot, 'DefaultFigureVisible', 'off');
        end
    end

    methods (Static, Access = private)
        function restore(previousPath, previousFolder, previousRng, ...
                previousVisibility, previousFigures, hadBackend, previousBackend)
            newFigures = setdiff(findall(groot, 'Type', 'figure'), previousFigures);
            delete(newFigures);
            set(groot, 'DefaultFigureVisible', previousVisibility);
            rng(previousRng);
            cd(previousFolder);
            path(previousPath);
            if hadBackend
                setappdata(groot, 'TIGREBackendState', previousBackend);
            elseif isappdata(groot, 'TIGREBackendState')
                rmappdata(groot, 'TIGREBackendState');
            end
        end
    end
end
