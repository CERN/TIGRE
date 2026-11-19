classdef BackendInitializationTest < backendtest.TestCase
    %BackendInitializationTest Public backend-selection contract tests.

    properties (TestParameter)
        Backend = struct('cpu', "cpu", 'metal', "metal")
    end

    methods (TestClassSetup)
        function addTigreToPath(testCase)
            matlabRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(matlabRoot));
        end
    end

    methods (Test)
        function testExplicitSelection(testCase, Backend)
            InitTIGRE('backend', Backend);

            actual = tigreBackend();

            testCase.verifyEqual(actual, Backend);
        end

        function testLogicalDevice(testCase, Backend)
            InitTIGRE('backend', Backend);

            gpuIds = GpuIds();

            testCase.verifyEqual(gpuIds.devices, int32(0));
            testCase.verifyEqual(length(gpuIds), 1);
        end

        function testDeviceNameCanBeSelected(testCase, Backend)
            InitTIGRE('backend', Backend);
            names = getGpuNames();

            gpuIds = GpuIds(names{1});

            testCase.verifyEqual(gpuIds.devices, int32(0));
        end

        function testSelectionSurvivesBaseWorkspaceClear(testCase)
            result = backendtest.runIsolated('clear', 'metal');

            testCase.verifyEqual(string(result.backend), "metal");
            testCase.verifyTrue(result.passed);
        end

        function testQueryBeforeInitializationErrors(testCase)
            rmappdata(groot, 'TIGREBackendState');

            action = @() tigreBackend();

            testCase.verifyError(action, 'TIGRE:Backend:NotInitialized');
        end

        function testInvalidBackendErrors(testCase)
            action = @() InitTIGRE('backend', 'invalid');

            testCase.verifyError( ...
                action, 'MATLAB:InputParser:ArgumentFailedValidation');
        end

        function testDefaultCudaFailsFastWhenUnavailable(testCase)
            testCase.assumeNotEqual(exist('Ax_mex', 'file'), 3);

            action = @() InitTIGRE();

            testCase.verifyError(action, 'TIGRE:InitTIGRE:BackendNotBuilt');
        end
    end
end
