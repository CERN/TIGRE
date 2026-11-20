# TIGRE website and documentation

Sources for https://cern.github.io/TIGRE/, built with Sphinx and
[pydata-sphinx-theme](https://pydata-sphinx-theme.readthedocs.io/).
Pages are written in Markdown (MyST) or reStructuredText.

## Build locally

```sh
pip install -r docs/requirements.txt
sphinx-build -W --keep-going -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html`. Building the site does not need
CUDA or a GPU, and does not need TIGRE to be installed.

## Publishing

- Every pull request that touches `docs/` builds the site with warnings
  treated as errors (`.github/workflows/docs.yml`). The built site is
  attached to the run as an artifact for preview.
- Every push to `master` deploys the site to GitHub Pages.

One-time repository setting: **Settings → Pages → Source: GitHub Actions**.

## Demos

The demo pages (`/demos/`) are generated from `Python/demos` and `MATLAB/Demos` by
[sphinx-gallery](https://sphinx-gallery.github.io/). The demo files are not
edited for the website: `docs/_ext/tigre_demos.py` writes converted copies to
`docs/_generated/demos_src` at build time, using the titles, descriptions,
sections and Python/MATLAB pairs listed in its `CATALOGUE`. **When you add or
rename a demo, add it to the catalogue.**

Running the demos needs TIGRE and a GPU, which the GitHub build does not have,
so outputs are rendered separately:

```sh
# on a GPU machine, from an up-to-date checkout of master
pip install -r docs/requirements.txt
python docs/tools/build_demos.py        # add --fresh to re-run every demo
git push --force origin docs-demos
# then Actions > Docs > Run workflow (master) to redeploy
```

The `docs-demos` branch holds only the rendered `docs/demos` folder (one
commit, replaced on every render). The Docs workflow unpacks it before
building, and sphinx-gallery reuses a rendered page as long as its converted
demo is unchanged. A demo changed after the last render shows code only until
the next render. Demos that need data or other packages (`NOT_RUN` in the
catalogue) and MATLAB-only demos are always shown as code.

## Gallery

`docs/gallery/index.md` is a hand-written page of reconstructions of real
scans (the old `Frontispiece/Gallery.md`). Most of its images are linked
rather than stored in the repository.

## API reference

- **Python** (`docs/api/python/`): `autodoc` + `napoleon` pages generated from
  the numpydoc docstrings in `Python/tigre`. The package is imported from the
  repository without being built; its CUDA extensions (and `h5py`, `tqdm`,
  `torch`) are mocked (`autodoc_mock_imports` in `conf.py`). The iterative
  algorithms' docstrings are assembled from shared pieces in
  `Python/tigre/algorithms/_docs.py`: edit them there.
- **MATLAB** (`docs/api/matlab/`): `docs/_ext/tigre_matlab_api.py` copies each
  function's help text (the first comment block, without the license) into
  `_generated/matlab/<group>.md`. The functions shown are listed in its
  `GROUPS`.
