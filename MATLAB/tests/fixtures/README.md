# CUDA reference fixture handoff

Two fixture sets are maintained:

- `cuda_reference.mat` covers the Stage 1 circular-geometry operators and
  TV family.
- `cuda_geometry_reference.mat` covers Stage 2 per-view geometry.

## Stage 2 one-command handoff

On the Windows machine with an NVIDIA GPU, check out the synchronized
implementation branch, start MATLAB in the repository root, and run:

```matlab
addpath(fullfile(pwd, 'MATLAB', 'tests', 'fixtures')); runCudaStage2Validation
```

The command builds and initializes CUDA, records MATLAB/computer/GPU/source
metadata, generates and validates the geometry fixture, runs `d13`, `d14`,
`d17`, `d18`, `d19`, and `d22` unchanged with figures hidden, and writes a
timestamped persistent log.

Each demo runs in a separate MATLAB process so its `clear` and `close all`
commands cannot affect the calling session. The validator restores the
caller's path, folder, RNG, backend selection, and figure visibility on
success or failure, and preserves existing figures. Failed
child runs retain their diagnostic files in the reported temporary folder.

Return both:

- `MATLAB/tests/fixtures/cuda_geometry_reference.mat`
- `MATLAB/tests/fixtures/cuda_stage2_validation_<timestamp>.log`

Do not regenerate or replace `cuda_reference.mat` during the Stage 2
handoff.

## Stage 1 fixture generation

To regenerate the original fixture:

1. Check out the implementation branch.
2. Start MATLAB in the repository root.
3. Run:

   ```matlab
   addpath(fullfile(pwd, 'MATLAB'))
   Compile('backend', 'cuda')
   addpath(fullfile(pwd, 'MATLAB', 'tests', 'fixtures'))
   generateCudaReferenceFixtures
   ```

4. Run the existing CUDA tests, then return:
   - `MATLAB/tests/fixtures/cuda_reference.mat`
   - the MATLAB command-window build/test log
   - the output of `version`, `computer`, and `getGpuNames()`

The fixtures are intentionally small and contain no source data. After
adding them to the branch, run `CudaReferenceTest` and
`CudaGeometryReferenceTest` on the Apple machine. Native-operation relative
L2 error must be at most `6e-3`.
