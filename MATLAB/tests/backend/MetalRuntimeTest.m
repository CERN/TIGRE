classdef MetalRuntimeTest < backendtest.TestCase
    %MetalRuntimeTest Initialisation can recover without unloading the MEX.

    properties (TestParameter)
        FailureMode = struct('missing', 'missing', 'corrupt', 'corrupt')
    end

    methods (Test)
        function testLibraryRecovery(testCase, FailureMode)
            result = backendtest.runIsolated('recovery', 'metal', FailureMode);

            testCase.verifyTrue(result.passed);
        end
    end
end
