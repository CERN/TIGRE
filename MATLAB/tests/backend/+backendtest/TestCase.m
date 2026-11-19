classdef TestCase < matlab.unittest.TestCase
    %TestCase Give each backend test a complete, reversible initialisation.

    methods (TestMethodSetup)
        function initialiseSession(testCase)
            testCase.applyFixture(backendtest.SessionFixture);
        end
    end
end
