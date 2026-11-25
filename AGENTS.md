# Instructions for AI coding agents

These notes are for AI coding tools (Claude Code, Codex, Copilot, Cursor and
others) working on TIGRE, and for the people using them. The full guide for all
contributors is [`contributing.md`](contributing.md); this file is the short
version, plus what agents most often get wrong here.

## What TIGRE is

A GPU toolbox for tomographic (CT) reconstruction, with the same features in
MATLAB and Python on top of shared CUDA code.

- `Common/CUDA/`: the GPU kernels (projectors, backprojectors, TV), shared by
  both languages. A change here affects MATLAB and Python.
- `Python/tigre/`: the Python package (`algorithms/`, `utilities/`, Cython
  bindings in `utilities/cuda_interface/`). Tests in `Python/tests/`, demos in
  `Python/demos/`.
- `MATLAB/`: the MATLAB toolbox (`Algorithms/`, `Utilities/`, `Demos/`,
  `Tests/`). `Compile.m` builds the MEX files.
- `docs/`: the website (Sphinx, MyST Markdown), https://cern.github.io/TIGRE/.
- `changelog.md`, `CITATION.cff`, `pyproject.toml`: release notes, citation and
  version.

## Building and testing

- Building needs an NVIDIA GPU, the CUDA toolkit and a C++ compiler:
  `pip install --group dev .` (re-run after changing CUDA or Cython code).
  `TIGRE_CUDA_ARCH=native` makes the build much faster while developing.
- `pytest` runs the core tests. Tests that need a GPU are marked
  `@pytest.mark.gpu` and are skipped without one; slow stress tests are marked
  `extended` (`pytest -m extended`).
- Most agents run without a GPU. Then:
  - write logic that can be tested on the CPU, and test it there (see the
    existing CPU-only tests, which load a module by file path or mock `Ax`/`Atb`);
  - say clearly in the pull request what you could **not** run (GPU tests,
    MATLAB), so a maintainer can run it.
- Docs: `pip install -r docs/requirements.txt`, then
  `sphinx-build -W --keep-going -b html docs docs/_build/html`. Warnings are
  errors in CI.

## Rules for changes

- **One thing per pull request.** Small, focused PRs get reviewed and merged;
  large mixed ones do not.
- **MATLAB and Python parity.** A feature added to one language goes into the
  other in the same PR, or the PR links an issue for it. Keep argument names
  and defaults the same in both.
- **Compatibility.** Do not rename, remove or change the meaning of anything
  public (functions, arguments, `geo` fields and their conventions) without the
  deprecation steps in the
  [compatibility policy](https://cern.github.io/TIGRE/user-guide/compatibility.html).
- **Changelog.** Add a line under `## Unreleased` in `changelog.md` for every
  change users can notice. Write it for users: what changed and what they will
  see, not how the code changed.
- **Tests: few, small, about behaviour.** Every test is compiled for and run
  on every PR, so each one has to earn its place.
  - Test what users rely on (a result, an error message, an argument), not
    internal details of how the code does it.
  - Before adding a test file, extend an existing test or parametrisation. A
    new file needs a reason.
  - Keep sizes small: volumes up to 64^3, a few iterations, a few seconds per
    test on a GPU. Anything bigger or rare is marked `extended`.
  - Not every bug fix needs a test. Add one when the bug could plausibly come
    back or affects results; then it should fail before the fix and pass after.
  - Say in the pull request why each new test is needed.
- **Line endings.** Many files use CRLF (most `.m` files, `setup.py`,
  `.gitignore`, parts of `README.md`). Keep each file's existing line endings;
  do not convert whole files.
- **Conventions.** Python arrays are indexed `(z, y, x)` for images and
  `(angle, v, u)` for projections; MATLAB projections are `proj(v, u, angle)`.
  At angle 0 the source is on +x. Angles are in radians; distances in any one
  unit, used consistently (the defaults and demos use mm).
- **Published methods.** If you implement a published algorithm, add its
  reference to both `TIGRE.bib` files (`Python/tigre/utilities/` and
  `MATLAB/Utilities/`) and cite it in the docstring.
- **Data.** Do not commit datasets or large binary files. Images are ignored by
  `.gitignore` except under `docs/**/images/` and `docs/_static/`.

## Pull requests

- Fill in the pull request template, including the test environment. If an
  agent wrote most of the change, say so, and say what a human checked.
- Do not push to `master`; work on a branch named after the change
  (`fix/fdk-offset-weights`, `docs/install-windows`).
- Keep the description short and factual: what changed, why, how it was
  tested, and anything you could not test.
