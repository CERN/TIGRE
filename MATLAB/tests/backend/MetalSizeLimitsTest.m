classdef MetalSizeLimitsTest < backendtest.TestCase
    %MetalSizeLimitsTest Reject unsafe output shapes without allocating them.

    properties (TestParameter)
        Projector = struct('siddon', 'Siddon', 'interpolated', 'interpolated')
        Backprojector = struct('fdk', 'FDK', 'matched', 'matched')
        Detector = struct('elementOverflow', [65536; 65536], ...
            'signedCoordinateOverflow', [2^31; 2], ...
            'hugeProduct', [2^31-1; 2^31-1])
        Volume = struct('elementOverflow', [2048; 2048; 1024], ...
            'signedCoordinateOverflow', [2^31; 1; 1])
    end

    methods (Test)
        function testProjectionSizeLimit(testCase, Projector, Detector)
            InitTIGRE('backend', 'metal');
            geometry = MetalSizeLimitsTest.geometry([2; 2; 2], Detector);

            testCase.verifyError(@() Ax(ones(2, 2, 2, 'single'), ...
                geometry, 0, Projector), 'TIGRE:Metal:SizeLimitExceeded');
        end

        function testBackprojectionSizeLimit(testCase, Backprojector, Volume)
            InitTIGRE('backend', 'metal');
            geometry = MetalSizeLimitsTest.geometry(Volume, [2; 2]);

            testCase.verifyError(@() Atb(ones(2, 2, 'single'), ...
                geometry, 0, Backprojector), 'TIGRE:Metal:SizeLimitExceeded');
        end

        function testProjectionCountIncludesAllViews(testCase)
            InitTIGRE('backend', 'metal');
            geometry = MetalSizeLimitsTest.geometry([2; 2; 2], [32768; 32768]);

            testCase.verifyError(@() Ax(ones(2, 2, 2, 'single'), ...
                geometry, [0, .1, .2, .3]), 'TIGRE:Metal:SizeLimitExceeded');
        end

        function testValidCallAfterSizeError(testCase)
            InitTIGRE('backend', 'metal');
            oversized = MetalSizeLimitsTest.geometry([2; 2; 2], [65536; 65536]);
            valid = MetalSizeLimitsTest.geometry([2; 2; 2], [3; 3]);
            testCase.verifyError(@() Ax(ones(2, 2, 2, 'single'), oversized, 0), ...
                'TIGRE:Metal:SizeLimitExceeded');

            actual = Ax(ones(2, 2, 2, 'single'), valid, 0);

            testCase.verifySize(actual, [3, 3]);
            testCase.verifyTrue(all(isfinite(actual), 'all'));
            testCase.verifyGreaterThan(norm(actual(:)), 0);
        end
    end

    methods (Static, Access = private)
        function value = geometry(volume, detector)
            value = struct('nVoxel', volume, 'sVoxel', volume, ...
                'dVoxel', ones(3, 1), 'nDetector', detector, ...
                'sDetector', detector, 'dDetector', ones(2, 1), ...
                'DSO', 40, 'DSD', 60);
        end
    end
end
