classdef BackendReconstructionTest < backendtest.TestCase
    %BackendReconstructionTest Quality and agreement on recorded CUDA data.

    properties (TestParameter)
        Backend = struct('cpu', "cpu", 'metal', "metal")
        GeometryCase = struct( ...
            'helical', struct('name', 'helical', 'errorLimit', .45), ...
            'originOffsetPerView', struct('name', 'originOffsetPerView', 'errorLimit', .45), ...
            'detectorRoll', struct('name', 'detectorRoll', 'errorLimit', .45), ...
            'arbitraryAngles', struct('name', 'arbitraryAngles', 'errorLimit', .45), ...
            'staticDetectorDbt', struct('name', 'staticDetectorDbt', 'errorLimit', .75))
        Workflow = struct( ...
            'helicalOSSART', struct('geometry', 'helical', 'algorithm', 'OS_SART'), ...
            'offsetSART', struct('geometry', 'originOffsetPerView', 'algorithm', 'SART'), ...
            'rotationOSSART', struct('geometry', 'detectorRoll', 'algorithm', 'OS_SART'), ...
            'arbitrarySIRT', struct('geometry', 'arbitraryAngles', 'algorithm', 'SIRT'), ...
            'dbtFDK', struct('geometry', 'staticDetectorDbt', 'algorithm', 'FDK'), ...
            'dbtSART', struct('geometry', 'staticDetectorDbt', 'algorithm', 'SART'))
    end

    methods (Test, TestTags = {'Reconstruction'})
        function testCglsQuality(testCase, Backend, GeometryCase)
            data = BackendReconstructionTest.fixtureCase(GeometryCase.name);
            InitTIGRE('backend', Backend);

            actual = CGLS(data.cuda.siddon, data.geometry, data.angles, ...
                10, 'Verbose', false);
            projected = Ax(actual, data.geometry, data.angles, 'Siddon');

            testCase.verifyClass(actual, 'single');
            testCase.verifySize(actual, size(data.image));
            testCase.verifyTrue(all(isfinite(actual), 'all'));
            testCase.verifyLessThanOrEqual(BackendReconstructionTest.relativeError( ...
                data.cuda.siddon, projected), .05, 'Projection residual');
            testCase.verifyLessThanOrEqual(BackendReconstructionTest.relativeError( ...
                data.image, actual), GeometryCase.errorLimit, 'Phantom recovery error');
        end

        function testCglsBackendAgreement(testCase, GeometryCase)
            data = BackendReconstructionTest.fixtureCase(GeometryCase.name);
            InitTIGRE('backend', 'cpu');
            reference = CGLS(data.cuda.siddon, data.geometry, data.angles, ...
                10, 'Verbose', false);
            InitTIGRE('backend', 'metal');

            actual = CGLS(data.cuda.siddon, data.geometry, data.angles, ...
                10, 'Verbose', false);

            testCase.verifyLessThanOrEqual( ...
                BackendReconstructionTest.relativeError(reference, actual), 1e-2);
        end

        function testWorkflowBackendAgreement(testCase, Workflow)
            data = BackendReconstructionTest.fixtureCase(Workflow.geometry);
            InitTIGRE('backend', 'cpu');
            reference = BackendReconstructionTest.reconstruct(data, Workflow.algorithm);
            InitTIGRE('backend', 'metal');

            actual = BackendReconstructionTest.reconstruct(data, Workflow.algorithm);

            testCase.verifyClass(actual, 'single');
            testCase.verifySize(actual, size(data.image));
            testCase.verifyTrue(all(isfinite(actual), 'all'));
            testCase.verifyGreaterThan(norm(double(reference(:))), 0);
            testCase.verifyLessThanOrEqual( ...
                BackendReconstructionTest.relativeError(reference, actual), 1e-2);
        end
    end

    methods (Static, Access = private)
        function value = fixtureCase(name)
            fixtureFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'fixtures', 'cuda_geometry_reference.mat');
            loaded = load(fixtureFile, 'fixture');
            validateCudaGeometryReferenceFixture(loaded.fixture);
            value = loaded.fixture.cases.(name);
        end

        function value = reconstruct(data, algorithm)
            rng(42, 'twister');
            switch algorithm
                case 'FDK'
                    value = FDK(data.cuda.interpolated, data.geometry, data.angles, ...
                        'parker', false, 'wang', false);
                case 'OS_SART'
                    value = OS_SART(data.cuda.interpolated, data.geometry, data.angles, ...
                        3, 'Verbose', false, 'BlockSize', 2, 'redundancy_weighting', false);
                case 'SART'
                    value = SART(data.cuda.interpolated, data.geometry, data.angles, ...
                        3, 'Verbose', false, 'OrderStrategy', 'ordered');
                case 'SIRT'
                    value = SIRT(data.cuda.interpolated, data.geometry, data.angles, ...
                        3, 'Verbose', false);
            end
        end

        function value = relativeError(reference, actual)
            value = norm(double(reference(:)) - double(actual(:))) / ...
                max(norm(double(reference(:))), eps);
        end
    end
end
