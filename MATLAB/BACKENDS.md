# MATLAB compute backends

TIGRE's MATLAB interface supports three session-wide native backends:

```matlab
InitTIGRE()                       % CUDA (legacy default)
InitTIGRE('backend', 'cuda')
InitTIGRE('backend', 'cpu')
InitTIGRE('backend', 'metal')
tigreBackend()                    % inspect the current selection
```

Calling `InitTIGRE` again switches the backend. The selection survives
`clear`, which allows the existing demos to run unchanged after
initialization. CPU and Metal expose one logical `GpuIds` device with ID 0.

Build a backend from the `MATLAB` folder or from any folder after adding it
to the MATLAB path:

```matlab
Compile('backend', 'cpu')
Compile('backend', 'metal')
Compile('backend', 'cuda')        % equivalent to Compile()
```

The CPU backend is portable C++17 and uses standard-library threads. The
Metal backend currently targets Apple Silicon macOS and builds both an
Objective-C++ MEX file and a Metal shader library. Generated binaries and
shader artifacts are intentionally not tracked.

## Geometry support

CPU and Metal currently support:

- cone/fan and parallel geometry;
- circular, helical, and arbitrary three-axis scan rotations;
- scalar or per-view DSO and DSD;
- scalar or per-view origin and detector offsets;
- scalar or per-view detector roll, pitch, and yaw;
- scalar or per-view centre-of-rotation correction;
- interpolated and Siddon forward projection;
- FDK, matched cone, and parallel backprojection;
- `minimizeTV`, `minimizeAwTV`, `im3DDenoise`, and `minimizeTV_prior`.

For the supported operations, CUDA output conventions are the compatibility
contract. CPU and Metal preserve CUDA detector-row ordering and its legacy
backprojection sampling and weighting behavior, including static-detector DBT
and combined parallel geometries.

Whole arrays must fit in host memory (CPU) or shared Metal memory; chunking
and multiple devices are outside this implementation. The CUDA-only pCT
projector is not available on CPU or Metal.

Metal also requires each array dimension to be below `2^31`, each volume or
complete projection stack to contain at most `2^32 - 1` elements, and each
buffer to fit within the GPU's `maxBufferLength`. These are indexing and
per-buffer limits, not a guarantee that enough memory is available. Unsafe
sizes are rejected before output allocation with
`TIGRE:Metal:SizeLimitExceeded`; reduce the problem size or use CPU/CUDA.
Buffer allocation can still fail if the device lacks available memory.

If the Metal shader library is missing or corrupt, the call reports
`TIGRE:Metal:LibraryNotFound`. Restore `tigre_kernels.metallib` beside the MEX
file (or rebuild Metal), then retry. A failed initialisation no longer leaves
a partially initialised runtime: restarting MATLAB or clearing the MEX is
not required after restoring the library.

On Apple Silicon, the Metal backend wraps page-aligned MATLAB input and output
storage directly as shared Metal buffers when both the address and byte length
meet Metal's no-copy requirements. It falls back to allocated shared buffers
and explicit copies for all other arrays. Calls remain synchronous, so wrapped
MATLAB storage remains valid until GPU execution completes.

The Metal Siddon kernel maintains the linear voxel address incrementally during
ray traversal. Compile `tigre_metal_mex.mm` with
`-DTIGRE_METAL_TIMING` in `CXXFLAGS` to print opt-in per-call timings for
geometry preparation, shared input setup, output allocation, host-observed GPU
execution, Metal GPU timestamps, and output copying. Standard builds contain
no timing calls or output.

## Validation

After building CPU and Metal, run validation from the repository root:

```matlab
testFolder = fullfile(pwd, 'MATLAB', 'tests', 'backend');
suite = testsuite(testFolder);
fastSuite = suite.selectIf(~matlab.unittest.selectors.HasTag('FullDemo'));
results = run(fastSuite);
assertSuccess(results);

fullDemoResults = runtests(testFolder, 'Tag', 'FullDemo');
assertSuccess(fullDemoResults);
```

Validation on Apple Silicon with MATLAB R2026a passed all 228 fast tests
and all 12 full demo/backend combinations after rebuilding CPU and Metal.
The combined fast run also preserved the caller's workspace sentinel,
existing figure, path, folder, RNG state, backend state, and figure visibility.
The standalone size-boundary tests passed with undefined-behaviour and
unsigned-overflow sanitizers enabled. The CUDA validation helper was checked
statically; running its CUDA build and fixture generation still requires an
NVIDIA machine.

Each fast test initialises TIGRE independently and restores the caller's
path, folder, random-number state, backend selection, and figure visibility.
Existing figures are preserved. Checks that deliberately clear workspaces,
simulate runtime failures, or run the six unchanged demos (`d13`, `d14`,
`d17`, `d18`, `d19`, `d22`) execute in separate `matlab -batch` processes.
The child runner requires the JVM and permission to start MATLAB processes;
it uses the same MATLAB installation as the parent. Full CPU demos can take
tens of minutes. A failed child retains its request, result (when available),
and log in the temporary folder printed by the test.

The fast suite covers session isolation, missing/corrupt shader recovery,
size rejection, operators, and reconstruction quality. Reconstruction tests
use the existing CUDA projections as inputs: CGLS must reduce relative
projection residuals to at most 5%, recover the phantom within geometry-specific
bounds (45%, or 75% for limited-angle DBT), and agree between CPU and Metal
within 1% relative L2 error. SART, OS-SART, SIRT, and FDK workflows also have
CPU/Metal agreement checks. The existing CUDA operator tolerance remains
`6e-3`; no NVIDIA hardware or additional test toolbox is needed to consume
the recorded fixtures.

Arithmetic boundaries can also be tested without MATLAB or a GPU:

```sh
c++ -std=c++17 -Wall -Wextra -Werror Common/Metal/tests/metal_limits_test.cpp -o /tmp/tigre_metal_limits_test
/tmp/tigre_metal_limits_test
```

The CUDA fixture handoff is documented in
`MATLAB/tests/fixtures/README.md`. Informational CPU/Metal timings are
available from `benchmarkBackends`. Reproducible Metal/CUDA `Ax` measurements,
including d13 and plot-free d13 CGLS, are available from
`benchmarkAxBackend`.
