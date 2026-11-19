function validateCudaGeometryReferenceFixture(fixture)
%VALIDATECUDAGEOMETRYREFERENCEFIXTURE Validate generated CUDA arrays.

requiredOutputs = {'interpolated', 'siddon', 'fdk', 'matched'};
caseNames = fieldnames(fixture.cases);
assert(~isempty(caseNames), ...
    'TIGRE:Tests:EmptyCudaGeometryFixture', ...
    'The CUDA geometry fixture has no cases.');

for caseIndex = 1:numel(caseNames)
    caseData = fixture.cases.(caseNames{caseIndex});
    validateArray(caseData.image, caseNames{caseIndex}, 'image');
    for outputIndex = 1:numel(requiredOutputs)
        outputName = requiredOutputs{outputIndex};
        assert(isfield(caseData.cuda, outputName), ...
            'TIGRE:Tests:MissingCudaGeometryOutput', ...
            'Case %s has no %s CUDA output.', ...
            caseNames{caseIndex}, outputName);
        validateArray( ...
            caseData.cuda.(outputName), caseNames{caseIndex}, outputName);
    end
end
end

function validateArray(value, caseName, arrayName)
assert(isnumeric(value) && ~isempty(value) && all(isfinite(value(:))), ...
    'TIGRE:Tests:InvalidCudaGeometryArray', ...
    'Case %s has an empty, nonnumeric, or nonfinite %s array.', ...
    caseName, arrayName);
end
