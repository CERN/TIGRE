classdef BackendFullDemoTest < matlab.unittest.TestCase
    %BackendFullDemoTest Run unchanged demos in independent MATLAB processes.

    properties (TestParameter)
        Backend = struct('cpu', "cpu", 'metal', "metal")
        Demo = struct( ...
            'helical', 'd13_HelicalGeometry.m', ...
            'offsets', 'd14_Offsets.m', ...
            'rotation', 'd17_DetectorRotation.m', ...
            'arbitraryAxis', 'd18_ArbitraryAxisOfRotation.m', ...
            'dbt', 'd19_DBTexample.m', ...
            'devices', 'd22_ListGPUs.m')
    end

    methods (Test, TestTags = {'FullDemo'})
        function testDemo(testCase, Backend, Demo)
            result = backendtest.runIsolated('demo', Backend, Demo);

            testCase.verifyTrue(result.passed);
            testCase.verifyEqual(string(result.backend), Backend);
            testCase.verifyNotEmpty(result.arrays);
        end
    end
end
