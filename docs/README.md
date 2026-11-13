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
