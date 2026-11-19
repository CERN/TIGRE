classdef BackendOperatorsTest < backendtest.TestCase
    %BackendOperatorsTest Native operator and general geometry tests.

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
        function testConeOperators(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendOperatorsTest.coneData();

            interpolated = Ax(image, geometry, angles, 'interpolated');
            siddon = Ax(image, geometry, angles, 'Siddon');
            fdk = Atb(interpolated, geometry, angles, 'FDK');
            matched = Atb(interpolated, geometry, angles, 'matched');

            testCase.verifySize(interpolated, [10 12 4]);
            testCase.verifySize(siddon, [10 12 4]);
            testCase.verifySize(fdk, [10 9 7]);
            testCase.verifySize(matched, [10 9 7]);
            testCase.verifyClass(interpolated, 'single');
            testCase.verifyTrue(all(isfinite( ...
                [interpolated(:); siddon(:); fdk(:); matched(:)])));
        end

        function testFanOperators(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendOperatorsTest.fanData();

            projections = Ax(image, geometry, angles, 'interpolated');
            reconstruction = Atb(projections, geometry, angles, 'matched');

            testCase.verifySize(projections, [1 14 5]);
            testCase.verifySize(reconstruction, [12 10]);
            testCase.verifyClass(reconstruction, 'single');
            testCase.verifyTrue(all(isfinite(reconstruction), 'all'));
        end

        function testParallelOperators(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendOperatorsTest.parallelData();

            interpolated = Ax(image, geometry, angles, 'interpolated');
            siddon = Ax(image, geometry, angles, 'Siddon');
            reconstruction = Atb(interpolated, geometry, angles, 'FDK');

            testCase.verifySize(interpolated, [7 13 4]);
            testCase.verifySize(siddon, [7 13 4]);
            testCase.verifySize(reconstruction, [11 9 7]);
            testCase.verifyTrue(all(isfinite( ...
                [interpolated(:); siddon(:); reconstruction(:)])));
        end

        function testTvFamily(testCase, Backend)
            InitTIGRE('backend', Backend);
            image = BackendOperatorsTest.tvImage();
            gpuIds = GpuIds();

            tv = minimizeTV(image, 0.05, 2, 'gpuids', gpuIds);
            adaptive = minimizeAwTV( ...
                image, 0.05, 2, -0.01, 'gpuids', gpuIds);
            denoised = im3DDenoise(image, 'TV', 2, 15, 'gpuids', gpuIds);
            piccs = minimizeTV_prior( ...
                image, 0.8 * image, 0.05, 2, 0.5, 'gpuids', gpuIds);

            testCase.verifySize(tv, size(image));
            testCase.verifySize(adaptive, size(image));
            testCase.verifySize(denoised, size(image));
            testCase.verifySize(piccs, size(image));
            testCase.verifyClass(tv, 'single');
            testCase.verifyTrue(all(isfinite( ...
                [tv(:); adaptive(:); denoised(:); piccs(:)])));
        end

        function testSupportsPerViewOffsets(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendOperatorsTest.coneData();
            geometry.offOrigin = [ ...
                sin(angles); cos(angles); linspace(-1, 1, numel(angles))];
            geometry.offDetector = [sin(angles); 0.5 * cos(angles)];
            geometry.COR = 0.25 * sin(angles);

            projections = Ax(image, geometry, angles, 'interpolated');

            testCase.verifySize(projections, [10 12 4]);
            testCase.verifyTrue(all(isfinite(projections), 'all'));
        end

        function testSupportsVaryingDistances(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendOperatorsTest.coneData();
            geometry.DSO = 40 + sin(angles);
            geometry.DSD = 60 + cos(angles);

            projections = Ax(image, geometry, angles, 'Siddon');

            testCase.verifySize(projections, [10 12 4]);
            testCase.verifyTrue(all(isfinite(projections), 'all'));
        end

        function testSupportsArbitraryRotation(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendOperatorsTest.coneData();
            angles = [angles; zeros(2, size(angles, 2))];
            angles(2, 2) = 0.1;
            geometry.rotDetector = [0.1; 0.05; -0.03];
            geometry.COR = 0.25;

            projections = Ax(image, geometry, angles, 'interpolated');

            testCase.verifySize(projections, [10 12 4]);
            testCase.verifyTrue(all(isfinite(projections), 'all'));
        end

        function testRejectsNonzeroLogicalDevice(testCase, Backend)
            InitTIGRE('backend', Backend);
            [image, geometry, angles] = BackendOperatorsTest.coneData();
            gpuIds = GpuIds();
            gpuIds.devices = int32(1);

            action = @() Ax( ...
                image, geometry, angles, 'interpolated', 'gpuids', gpuIds);

            testCase.verifyError(action, 'TIGRE:Backend:InvalidDevice');
        end
    end

    methods (Test)
        function testInterpolatedCpuMetalAgreement(testCase)
            [image, geometry, angles] = BackendOperatorsTest.coneData();

            cpu = BackendOperatorsTest.runAx( ...
                "cpu", image, geometry, angles, 'interpolated');
            metal = BackendOperatorsTest.runAx( ...
                "metal", image, geometry, angles, 'interpolated');

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testSiddonCpuMetalAgreement(testCase)
            [image, geometry, angles] = BackendOperatorsTest.coneData();

            cpu = BackendOperatorsTest.runAx( ...
                "cpu", image, geometry, angles, 'Siddon');
            metal = BackendOperatorsTest.runAx( ...
                "metal", image, geometry, angles, 'Siddon');

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testFdkCpuMetalAgreement(testCase)
            [image, geometry, angles] = BackendOperatorsTest.coneData();
            InitTIGRE('backend', 'cpu');
            projections = Ax(image, geometry, angles, 'interpolated');

            cpu = BackendOperatorsTest.runAtb( ...
                "cpu", projections, geometry, angles, 'FDK');
            metal = BackendOperatorsTest.runAtb( ...
                "metal", projections, geometry, angles, 'FDK');

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testMatchedCpuMetalAgreement(testCase)
            [image, geometry, angles] = BackendOperatorsTest.coneData();
            InitTIGRE('backend', 'cpu');
            projections = Ax(image, geometry, angles, 'interpolated');

            cpu = BackendOperatorsTest.runAtb( ...
                "cpu", projections, geometry, angles, 'matched');
            metal = BackendOperatorsTest.runAtb( ...
                "metal", projections, geometry, angles, 'matched');

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testTvCpuMetalAgreement(testCase)
            image = BackendOperatorsTest.tvImage();

            cpu = BackendOperatorsTest.runTv("cpu", image);
            metal = BackendOperatorsTest.runTv("metal", image);

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testAdaptiveTvCpuMetalAgreement(testCase)
            image = BackendOperatorsTest.tvImage();

            cpu = BackendOperatorsTest.runAdaptiveTv("cpu", image);
            metal = BackendOperatorsTest.runAdaptiveTv("metal", image);

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testDenoiseCpuMetalAgreement(testCase)
            image = BackendOperatorsTest.tvImage();

            cpu = BackendOperatorsTest.runDenoise("cpu", image);
            metal = BackendOperatorsTest.runDenoise("metal", image);

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testPiccsCpuMetalAgreement(testCase)
            image = BackendOperatorsTest.tvImage();

            cpu = BackendOperatorsTest.runPiccs("cpu", image);
            metal = BackendOperatorsTest.runPiccs("metal", image);

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpu, metal), 5e-3);
        end

        function testGeneralConeGeometryCpuMetalAgreement(testCase)
            [image, geometry, angles] = ...
                BackendOperatorsTest.generalGeometryData('cone');

            cpuInterpolated = BackendOperatorsTest.runAx( ...
                "cpu", image, geometry, angles, 'interpolated');
            metalInterpolated = BackendOperatorsTest.runAx( ...
                "metal", image, geometry, angles, 'interpolated');
            cpuSiddon = BackendOperatorsTest.runAx( ...
                "cpu", image, geometry, angles, 'Siddon');
            metalSiddon = BackendOperatorsTest.runAx( ...
                "metal", image, geometry, angles, 'Siddon');
            cpuFdk = BackendOperatorsTest.runAtb( ...
                "cpu", cpuInterpolated, geometry, angles, 'FDK');
            metalFdk = BackendOperatorsTest.runAtb( ...
                "metal", cpuInterpolated, geometry, angles, 'FDK');
            cpuMatched = BackendOperatorsTest.runAtb( ...
                "cpu", cpuInterpolated, geometry, angles, 'matched');
            metalMatched = BackendOperatorsTest.runAtb( ...
                "metal", cpuInterpolated, geometry, angles, 'matched');

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError( ...
                    cpuInterpolated, metalInterpolated), 5e-3);
            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpuSiddon, metalSiddon), ...
                5e-3);
            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError(cpuFdk, metalFdk), 5e-3);
            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError( ...
                    cpuMatched, metalMatched), 5e-3);
        end

        function testGeneralParallelGeometryCpuMetalAgreement(testCase)
            [image, geometry, angles] = ...
                BackendOperatorsTest.generalGeometryData('parallel');

            cpuProjection = BackendOperatorsTest.runAx( ...
                "cpu", image, geometry, angles, 'interpolated');
            metalProjection = BackendOperatorsTest.runAx( ...
                "metal", image, geometry, angles, 'interpolated');
            cpuBackprojection = BackendOperatorsTest.runAtb( ...
                "cpu", cpuProjection, geometry, angles, 'FDK');
            metalBackprojection = BackendOperatorsTest.runAtb( ...
                "metal", cpuProjection, geometry, angles, 'FDK');

            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError( ...
                    cpuProjection, metalProjection), 5e-3);
            testCase.verifyLessThanOrEqual( ...
                BackendOperatorsTest.relativeError( ...
                    cpuBackprojection, metalBackprojection), 5e-3);
        end
    end

    methods (Static, Access = private)
        function [image, geometry, angles] = coneData()
            geometry = BackendOperatorsTest.geometry( ...
                [10; 9; 7], [12; 10], 'cone');
            angles = linspace(0, 2 * pi, 5);
            angles(end) = [];
            image = BackendOperatorsTest.phantom(geometry.nVoxel);
        end

        function [image, geometry, angles] = fanData()
            geometry = BackendOperatorsTest.geometry( ...
                [12; 10; 1], [14; 1], 'cone');
            angles = linspace(0, 2 * pi, 6);
            angles(end) = [];
            image = BackendOperatorsTest.phantom(geometry.nVoxel);
        end

        function [image, geometry, angles] = parallelData()
            geometry = BackendOperatorsTest.geometry( ...
                [11; 9; 7], [13; 7], 'parallel');
            angles = linspace(0, 2 * pi, 5);
            angles(end) = [];
            image = BackendOperatorsTest.phantom(geometry.nVoxel);
        end

        function [image, geometry, angles] = generalGeometryData(mode)
            geometry = BackendOperatorsTest.geometry( ...
                [12; 11; 9], [14; 11], mode);
            alpha = linspace(0, 2 * pi, 7);
            alpha(end) = [];
            angles = [ ...
                alpha; ...
                0.1 * sin(alpha); ...
                0.06 * cos(alpha)];
            geometry.offOrigin = [ ...
                sin(alpha); cos(alpha); linspace(-1.5, 1.5, numel(alpha))];
            geometry.offDetector = [ ...
                1.5 * sin(alpha); 0.75 * cos(alpha)];
            geometry.COR = 0.2 * sin(alpha);
            geometry.rotDetector = [ ...
                0.08 * sin(alpha); ...
                0.04 * cos(alpha); ...
                0.06 * sin(2 * alpha)];
            geometry.DSO = 40 + sin(alpha);
            geometry.DSD = 60 + cos(alpha);
            image = BackendOperatorsTest.phantom(geometry.nVoxel);
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
            geometry.accuracy = 0.5;
        end

        function image = phantom(dimensions)
            [x, y, z] = ndgrid( ...
                linspace(-1, 1, dimensions(1)), ...
                linspace(-1, 1, dimensions(2)), ...
                linspace(-1, 1, dimensions(3)));
            image = single((x.^2 + 1.2 * y.^2 + 0.8 * z.^2) < 0.55);
        end

        function image = tvImage()
            image = reshape(single(mod(0:(10 * 9 * 7 - 1), 23)) / 23, ...
                [10 9 7]);
        end

        function result = runAx(backend, image, geometry, angles, projector)
            InitTIGRE('backend', backend);
            result = Ax(image, geometry, angles, projector);
        end

        function result = runAtb(backend, projections, geometry, angles, projector)
            InitTIGRE('backend', backend);
            result = Atb(projections, geometry, angles, projector);
        end

        function result = runTv(backend, image)
            InitTIGRE('backend', backend);
            result = minimizeTV(image, 0.05, 2, 'gpuids', GpuIds());
        end

        function result = runAdaptiveTv(backend, image)
            InitTIGRE('backend', backend);
            result = minimizeAwTV( ...
                image, 0.05, 2, -0.01, 'gpuids', GpuIds());
        end

        function result = runDenoise(backend, image)
            InitTIGRE('backend', backend);
            result = im3DDenoise(image, 'TV', 2, 15, 'gpuids', GpuIds());
        end

        function result = runPiccs(backend, image)
            InitTIGRE('backend', backend);
            result = minimizeTV_prior( ...
                image, 0.8 * image, 0.05, 2, 0.5, 'gpuids', GpuIds());
        end

        function value = relativeError(reference, actual)
            value = norm(double(reference(:) - actual(:))) / ...
                max(norm(double(reference(:))), eps);
        end
    end
end
