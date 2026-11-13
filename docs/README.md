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

## Gallery

The gallery is generated from `Python/demos` and `MATLAB/Demos` by
[sphinx-gallery](https://sphinx-gallery.github.io/). The demo files are not
edited for the website: `docs/_ext/tigre_gallery.py` writes converted copies to
`docs/_generated/gallery_src` at build time, using the titles, descriptions,
sections and Python/MATLAB pairs listed in its `CATALOGUE`. **When you add or
rename a demo, add it to the catalogue.**

Running the demos needs TIGRE and a GPU, which the GitHub build does not have,
so outputs are rendered separately:

```sh
# on a GPU machine, from an up-to-date checkout of master
pip install -r docs/requirements.txt
python docs/tools/build_gallery.py        # add --fresh to re-run every demo
git push --force origin docs-gallery
# then Actions > Docs > Run workflow (master) to redeploy
```

The `docs-gallery` branch holds only the rendered `docs/gallery` folder (one
commit, replaced on every render). The Docs workflow unpacks it before
building, and sphinx-gallery reuses a rendered page as long as its converted
demo is unchanged. A demo changed after the last render shows code only until
the next render. Demos that need data or other packages (`NOT_RUN` in the
catalogue) and MATLAB-only demos are always shown as code.
