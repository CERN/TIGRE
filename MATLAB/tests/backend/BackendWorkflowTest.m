classdef BackendWorkflowTest < backendtest.TestCase
    %BackendWorkflowTest Reduced public workflows based on TIGRE demos.

    properties (TestParameter)
        Backend = struct('cpu', "cpu", 'metal', "metal")
    end

    methods (TestClassSetup)
        function addTigreToPath(testCase)
            matlabRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(matlabRoot));
        end
    end

    methods (TestMethodSetup)
        function seedRandomStream(testCase)
            previous = rng;
            testCase.addTeardown(@() rng(previous));
            rng(42, 'twister');
        end
    end

    methods (Test)
        function testProjectionNoiseFdkWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendWorkflowTest.coneData();

            projections = Ax(image, geometry, angles, 'interpolated');
            noisy = addCTnoise( ...
                projections, 'Poisson', 1e4, 'Gaussian', [0 0.1]);
            reconstruction = FDK( ...
                noisy, geometry, angles, 'parker', false, 'wang', false);

            testCase.verifySize(reconstruction, size(image));
            testCase.verifyClass(reconstruction, 'single');
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testCglsWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            projections = Ax(image, geometry, angles, 'Siddon');

            reconstruction = CGLS( ...
                projections, geometry, angles, 1, 'verbose', false);

            testCase.verifySize(reconstruction, size(image));
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testParallelOsSartWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendWorkflowTest.parallelData();
            projections = Ax(image, geometry, angles, 'interpolated');

            reconstruction = OS_SART( ...
                projections, geometry, angles, 1, ...
                'BlockSize', 2, 'Verbose', false, ...
                'redundancy_weighting', false);

            testCase.verifySize(reconstruction, size(image));
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testAdaptiveTvWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            projections = Ax(image, geometry, angles, 'interpolated');

            reconstruction = AwASD_POCS( ...
                projections, geometry, angles, 1, ...
                'TViter', 1, 'alpha', 0.001, 'maxL2err', inf, ...
                'Verbose', false, 'redundancy_weighting', false);

            testCase.verifySize(reconstruction, size(image));
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testDenoisingWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            image = BackendWorkflowTest.tvImage();

            denoised = im3DDenoise(image, 'TV', 3, 15);

            testCase.verifySize(denoised, size(image));
            testCase.verifyTrue(all(isfinite(denoised), 'all'));
        end

        function testPiccsWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            projections = Ax(image, geometry, angles, 'interpolated');

            reconstruction = PICCS( ...
                projections, geometry, angles, 1, image, ...
                'TViter', 1, 'alpha', 0.001, 'maxL2err', inf, ...
                'Verbose', false);

            testCase.verifySize(reconstruction, size(image));
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testHelicalWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendWorkflowTest.helicalData();
            projections = Ax(image, geometry, angles, 'interpolated');

            osSart = OS_SART( ...
                projections, geometry, angles, 1, ...
                'BlockSize', 2, 'Verbose', false, ...
                'redundancy_weighting', false);
            cgls = CGLS( ...
                projections, geometry, angles, 1, 'Verbose', false);

            testCase.verifySize(osSart, size(image));
            testCase.verifySize(cgls, size(image));
            testCase.verifyTrue(all(isfinite([osSart(:); cgls(:)])));
        end

        function testOffsetWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendWorkflowTest.offsetData();
            projections = Ax(image, geometry, angles, 'interpolated');

            reconstruction = SART( ...
                projections, geometry, angles, 1, 'Verbose', false);

            testCase.verifySize(reconstruction, size(image));
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testRotatedDetectorWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = ...
                BackendWorkflowTest.rotatedDetectorData();
            projections = Ax(image, geometry, angles, 'interpolated');

            reconstruction = OS_SART( ...
                projections, geometry, angles, 1, ...
                'BlockSize', 2, 'Verbose', false, ...
                'redundancy_weighting', false);

            testCase.verifySize(reconstruction, size(image));
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testArbitraryAxisWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = ...
                BackendWorkflowTest.arbitraryAxisData();
            projections = Ax(image, geometry, angles, 'Siddon');

            sirt = SIRT( ...
                projections, geometry, angles, 1, 'Verbose', false);
            cgls = CGLS( ...
                projections, geometry, angles, 1, 'Verbose', false);

            testCase.verifySize(sirt, size(image));
            testCase.verifySize(cgls, size(image));
            testCase.verifyTrue(all(isfinite([sirt(:); cgls(:)])));
        end

        function testStaticDetectorWorkflow(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = ...
                BackendWorkflowTest.staticDetectorData();
            projections = Ax(image, geometry, angles, 'interpolated');

            fdk = FDK( ...
                projections, geometry, angles, ...
                'parker', false, 'wang', false);
            sart = SART( ...
                projections, geometry, angles, 1, ...
                'OrderStrategy', 'ordered', 'Verbose', false);

            testCase.verifySize(fdk, size(image));
            testCase.verifySize(sart, size(image));
            testCase.verifyTrue(all(isfinite([fdk(:); sart(:)])));
        end

        function testFdkCpuMetalAgreement(testCase)
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            InitTIGRE('backend', 'cpu');
            projections = Ax(image, geometry, angles, 'interpolated');

            cpu = BackendWorkflowTest.runFdk( ...
                "cpu", projections, geometry, angles);
            metal = BackendWorkflowTest.runFdk( ...
                "metal", projections, geometry, angles);

            testCase.verifyLessThanOrEqual( ...
                BackendWorkflowTest.nrmse(cpu, metal), 1e-2);
        end

        function testCglsCpuMetalAgreement(testCase)
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            InitTIGRE('backend', 'cpu');
            projections = Ax(image, geometry, angles, 'Siddon');

            cpu = BackendWorkflowTest.runCgls( ...
                "cpu", projections, geometry, angles);
            metal = BackendWorkflowTest.runCgls( ...
                "metal", projections, geometry, angles);

            testCase.verifyLessThanOrEqual( ...
                BackendWorkflowTest.nrmse(cpu, metal), 1e-2);
        end
    end

    methods (Static, Access = private)
        function [image, geometry, angles] = coneData()
            geometry = BackendWorkflowTest.geometry( ...
                [12; 12; 8], [16; 10], 'cone');
            angles = linspace(0, 2 * pi, 9);
            angles(end) = [];
            image = BackendWorkflowTest.phantom(geometry.nVoxel);
        end

        function [image, geometry, angles] = parallelData()
            geometry = BackendWorkflowTest.geometry( ...
                [12; 12; 8], [16; 8], 'parallel');
            angles = linspace(0, pi, 9);
            angles(end) = [];
            image = BackendWorkflowTest.phantom(geometry.nVoxel);
        end

        function [image, geometry, angles] = helicalData()
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            geometry.offOrigin = zeros(3, numel(angles));
            geometry.offOrigin(3, :) = linspace(-2, 2, numel(angles));
        end

        function [image, geometry, angles] = offsetData()
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            geometry.offOrigin = [ ...
                1.5 * sin(angles); ...
                linspace(-1, 1, numel(angles)); ...
                cos(angles)];
            geometry.offDetector = [ ...
                2 * sin(angles); 1.5 * cos(angles)];
        end

        function [image, geometry, angles] = rotatedDetectorData()
            [image, geometry, angles] = BackendWorkflowTest.coneData();
            geometry.rotDetector = [ ...
                0.1 * sin(angles); ...
                0.05 * cos(angles); ...
                0.08 * sin(2 * angles)];
        end

        function [image, geometry, angles] = arbitraryAxisData()
            [image, geometry, alpha] = BackendWorkflowTest.coneData();
            angles = [ ...
                alpha; ...
                0.12 * sin(alpha); ...
                0.08 * cos(alpha)];
        end

        function [image, geometry, angles] = staticDetectorData()
            geometry = BackendWorkflowTest.geometry( ...
                [10; 12; 8], [14; 10], 'cone');
            geometry.DSO = 50;
            geometry.DSD = 65;
            geometry.offOrigin = [2; 0; 1];
            geometry.offDetector = [0; 2];
            angles = linspace(-pi / 12, pi / 12, 7);
            geometry = staticDetectorGeo(geometry, angles);
            image = BackendWorkflowTest.phantom(geometry.nVoxel);
        end

        function geometry = geometry(nVoxel, nDetector, mode)
            geometry = defaultGeometry( ...
                'nVoxel', nVoxel, 'nDetector', nDetector, 'mode', mode);
            geometry.nVoxel = double(nVoxel);
            geometry.nDetector = double(nDetector);
            geometry.sVoxel = double(nVoxel);
            geometry.dVoxel = geometry.sVoxel ./ geometry.nVoxel;
            geometry.dDetector = [1; 1];
            geometry.sDetector = geometry.nDetector .* geometry.dDetector;
            geometry.DSO = 40;
            geometry.DSD = 60;
            geometry.accuracy = 0.75;
        end

        function image = phantom(dimensions)
            [x, y, z] = ndgrid( ...
                linspace(-1, 1, dimensions(1)), ...
                linspace(-1, 1, dimensions(2)), ...
                linspace(-1, 1, dimensions(3)));
            image = single((x.^2 + y.^2 + z.^2) < 0.55);
        end

        function image = tvImage()
            image = reshape(single(mod(0:(12 * 12 * 8 - 1), 29)) / 29, ...
                [12 12 8]);
        end

        function reconstruction = runFdk(backend, projections, geometry, angles)
            InitTIGRE('backend', backend);
            reconstruction = FDK( ...
                projections, geometry, angles, ...
                'parker', false, 'wang', false);
        end

        function reconstruction = runCgls(backend, projections, geometry, angles)
            InitTIGRE('backend', backend);
            reconstruction = CGLS( ...
                projections, geometry, angles, 1, 'verbose', false);
        end

        function value = nrmse(reference, actual)
            value = norm(double(reference(:) - actual(:))) / ...
                max(norm(double(reference(:))), eps);
        end
    end
end
