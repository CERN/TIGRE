# Contributing to TIGRE

Thanks for helping improve TIGRE! Bug reports, fixes, new algorithms, demos and
documentation are all welcome. This guide explains how to get set up and what we
look for in a pull request.

By participating you agree to follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## Asking questions and reporting bugs

- **Questions about using TIGRE** (geometry setup, which algorithm to use, loading
  your scanner data): please open a [Discussion](https://github.com/CERN/TIGRE/discussions)
  rather than an issue. Check the [FAQ](Frontispiece/FAQ.md) first.
- **Bugs**: open an [issue](https://github.com/CERN/TIGRE/issues/new/choose) using the
  bug template. Include a minimal script that reproduces the problem, and your OS,
  GPU, CUDA version and MATLAB/Python version. Most TIGRE bugs depend on these.
- **Feature ideas**: open an issue with the feature template so we can discuss the
  design before you write a lot of code.

## Setting up a development environment

You need an NVIDIA GPU, the CUDA toolkit and a C++ compiler (MSVC on Windows,
gcc on Linux). See the full installation guides for
[Python](Frontispiece/python_installation.md) and
[MATLAB](Frontispiece/MATLAB_installation.md).

### Python

```bash
git clone https://github.com/<your-username>/TIGRE.git
cd TIGRE
python -m venv .venv            # or a conda environment
source .venv/bin/activate        # .venv\Scripts\activate on Windows
pip install --group dev .        # builds the CUDA extensions and installs dev tools
```

Re-run `pip install .` after changing any CUDA (`Common/CUDA`) or Cython
(`Python/tigre/utilities/cuda_interface`) code, since those need recompiling.

### MATLAB

Run `MATLAB/Compile.m` once to build the MEX files, then `InitTIGRE` at the start of
each session to add TIGRE to the path.

## Running the tests

```bash
pytest
```

Tests live in `Python/tests/`. Tests that need a GPU are marked with
`@pytest.mark.gpu` (or `pytestmark = pytest.mark.gpu` for a whole file) and are
skipped automatically on machines without one, which is how they run on GitHub
Actions. Please run the full suite on a GPU machine before submitting changes to
projectors, kernels or algorithms.

MATLAB tests are in `MATLAB/Tests/`.

## Making a change

1. Fork the repository and create a branch from `master` with a descriptive name
   (`fix/fdk-offset-weights`, `feature/new-loader`).
2. Keep each pull request focused on one thing. Small PRs get reviewed faster.
3. Add tests for new behaviour and for bug fixes where possible. Pure-Python logic
   should be testable without a GPU; mock `Ax`/`Atb` if needed (see
   `Python/tests/test_krylov_restart_and_norms.py` for an example).
4. Update the demos and documentation if you change user-facing behaviour.
5. Add a line to [`changelog.md`](changelog.md) under the next release.
6. Open the pull request and fill in the template.

### Code style

- Python is formatted with [black](https://black.readthedocs.io/) (line length 100)
  and checked with flake8. Run `black .` before committing.
- Use numpy-style docstrings for public functions, and describe every parameter.
- MATLAB code should follow the style of the surrounding files, with a help block
  at the top of each function.

### MATLAB and Python parity

TIGRE aims to offer the same features in both languages. If you add a feature or
algorithm to one, either port it to the other in the same PR or open an issue so
someone else can. Keep argument names and defaults the same across both.

### Adding an algorithm

- **Python**: add it to the matching family module in `Python/tigre/algorithms/`
  (for example `krylov_subspace_algorithms.py`, `pocs_algorithms.py`), deriving from
  `IterativeReconAlg` where possible, and export it in `Python/tigre/algorithms/__init__.py`.
- **MATLAB**: add the function to `MATLAB/Algorithms/`.
- Add or extend a demo showing how to use it.
- If it implements a published method, add the reference to `TIGRE.bib`
  (both `Python/tigre/utilities/TIGRE.bib` and `MATLAB/Utilities/TIGRE.bib`) and
  cite it in the docstring.

### Changing CUDA kernels

The kernels in `Common/CUDA` are shared by both languages, so changes affect
MATLAB and Python. Please mention in the PR which GPU(s) and CUDA version(s) you
tested on, and include timings if the change is meant to improve performance.

## Recognition

We use [all-contributors](https://allcontributors.org/) to thank everyone who
contributes, including code, documentation, bug reports and user testing.
Ask in your PR to be added.

## License

By contributing to TIGRE, you agree that your contributions will be licensed
under its [BSD 3-Clause license](LICENSE.txt).
