classdef BackendIsolationTest < backendtest.TestCase
    %BackendIsolationTest Child checks and fixtures preserve the caller.

    methods (Test)
        function testSuccessfulChildPreservesSession(testCase)
            sentinelName = BackendIsolationTest.baseSentinel(testCase);
            existingFigure = figure('Visible', 'off');
            before = BackendIsolationTest.sessionState();

            result = backendtest.runIsolated('clear', 'metal');

            testCase.verifyTrue(result.passed);
            testCase.verifyEqual(evalin('base', sentinelName), 42);
            testCase.verifyTrue(isgraphics(existingFigure));
            testCase.verifyEqual(BackendIsolationTest.sessionState(), before);
        end

        function testFailedChildPreservesSession(testCase)
            sentinelName = BackendIsolationTest.baseSentinel(testCase);
            existingFigure = figure('Visible', 'off');
            before = BackendIsolationTest.sessionState();

            testCase.verifyError(@() backendtest.runIsolated('fail', 'cpu'), ...
                'TIGRE:Tests:ChildFailed');

            testCase.verifyEqual(evalin('base', sentinelName), 42);
            testCase.verifyTrue(isgraphics(existingFigure));
            testCase.verifyEqual(BackendIsolationTest.sessionState(), before);
        end

        function testMissingResultFails(testCase)
            testCase.verifyError(@() backendtest.runIsolated('noResult', 'cpu'), ...
                'TIGRE:Tests:MissingChildResult');
        end

        function testFixtureRestoresSessionAfterError(testCase)
            existingFigure = figure('Visible', 'off');
            before = BackendIsolationTest.sessionState();

            testCase.verifyError(@BackendIsolationTest.failingFixture, ...
                'MATLAB:unittest:Fixture:ApplyAndRunFailed');

            testCase.verifyTrue(isgraphics(existingFigure));
            testCase.verifyEqual(BackendIsolationTest.sessionState(), before);
        end

        function testFixtureRestoresUninitialisedSession(testCase)
            rmappdata(groot, 'TIGREBackendState');
            before = BackendIsolationTest.sessionState();

            testCase.verifyError(@BackendIsolationTest.failingFixture, ...
                'MATLAB:unittest:Fixture:ApplyAndRunFailed');

            testCase.verifyFalse(isappdata(groot, 'TIGREBackendState'));
            testCase.verifyEqual(BackendIsolationTest.sessionState(), before);
        end
    end

    methods (Static, Access = private)
        function name = baseSentinel(testCase)
            [~, suffix] = fileparts(tempname);
            name = matlab.lang.makeValidName(['tigreIsolation_', suffix]);
            assignin('base', name, 42);
            testCase.addTeardown(@() evalin('base', ['clear ', name]));
        end

        function state = sessionState()
            state = struct('path', path, 'folder', pwd, 'rng', rng, ...
                'visibility', get(groot, 'DefaultFigureVisible'), ...
                'figures', findall(groot, 'Type', 'figure'), ...
                'hasBackend', isappdata(groot, 'TIGREBackendState'), ...
                'backend', getappdata(groot, 'TIGREBackendState'));
        end

        function failingFixture()
            fixture = backendtest.SessionFixture;
            % Let the framework execute registered teardown callbacks before
            % rethrowing, rather than calling Fixture.teardown directly.
            fixture.applyAndRun(@BackendIsolationTest.mutateAndFail);
        end

        function mutateAndFail()
            InitTIGRE('backend', 'metal');
            rng(17, 'twister');
            cd(tempdir);
            set(groot, 'DefaultFigureVisible', 'on');
            figure('Visible', 'off');
            error('TIGRE:Tests:IntentionalFixtureFailure', 'Intentional fixture probe.');
        end
    end
end
