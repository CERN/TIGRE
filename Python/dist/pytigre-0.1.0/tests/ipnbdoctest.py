import os
import subprocess
import tempfile

import nbformat


def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """
    dirname, __ = os.path.split(path)
    os.chdir(dirname)
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=60",
                "--output", fout.name, path]
        subprocess.check_call(args)

        fout.seek(0)
        nb = nbformat.read(fout, nbformat.current_nbformat)

    errors = [output for cell in nb.cells if "outputs" in cell
              for output in cell["outputs"] \
              if output.output_type == "error"]

    return nb, errors


def test_ipynb():
    nb, errors = _notebook_run('../tigre/demos/d01_create_geometry.ipynb')
    assert errors == []


if __name__ == '__main__':
    test_ipynb()
