"""
End-to-end smoke tests that invoke tigre_reconstruction.py exactly as a user
would from the shell: a real subprocess, real argv parsing, the real fixture
tilt series in this folder, and a real GPU reconstruction via TIGRE (method
wbp only). This is the complement to test_tigre_reconstruction.py, which
calls the class directly with the heavy TIGRE/mrcfile/torch calls mocked out
and therefore never proves the command line itself works.

Fixtures used (all in this tests/ folder):
  ts.mrcs -- real tilt series (29 images, 256x256)
  ts.tlt  -- matching tilt angles (29 values)
  ts.xf   -- 29 identity (unitary) transforms: "1.0 0.0 0.0 1.0 0.0 0.0" per
             line, i.e. no rotation/scale/shift, generated once alongside
             ts.mrcs/ts.tlt so tigre_reconstruction.py's alignment step is a
             no-op.

--tigreInterpolation is intentionally never used: this file only tests the
default path where alignment is applied via applyTorchTransforms.

Requires a CUDA GPU (TIGRE's backprojection runs on it) and the `tigre`
conda env (torch, mrcfile, emtable, tigre all installed there).
"""
import os
import subprocess
import sys

import numpy as np
import mrcfile

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
GITTIGRE_ROOT = os.path.dirname(TESTS_DIR)
SCRIPT = os.path.join(GITTIGRE_ROOT, 'tigre_reconstruction.py')

TS_PATH = os.path.join(TESTS_DIR, 'ts.mrcs')
ANGLES_PATH = os.path.join(TESTS_DIR, 'ts.tlt')
XF_PATH = os.path.join(TESTS_DIR, 'ts.xf')
TOMO_THICKNES = 64

ITERATIONS_GRADIENT = 20
ITERATIONS_KRYLOV = 30
ITERATIONS_STATISTICAL = 20
ITERATIONS_VARIATIONAL = 20

GPU_ID = '0'

with mrcfile.open(TS_PATH, permissive=True, header_only=True) as _mrc:
    N_IMAGES, H, W = _mrc.header.nz.item(), _mrc.header.ny.item(), _mrc.header.nx.item()


def runCli(extraArgs, timeout=3600):
    return subprocess.run(
        [sys.executable, SCRIPT] + extraArgs,
        cwd=GITTIGRE_ROOT,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def assertReconstructionOk(result, outPath, expectedShape=(TOMO_THICKNES, H, W)):
    assert result.returncode == 0, result.stdout + result.stderr
    assert outPath.exists()
    with mrcfile.open(str(outPath), permissive=True) as mrc:
        data = mrc.data
        assert data.shape == expectedShape
        assert np.isfinite(data).all()


class TestWbp:

    def test_wbp_end_to_end(self, tmp_path):
        outPath = tmp_path / 'wbp.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'wbp', '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_wbp_custom_filter_end_to_end(self, tmp_path):
        outPath = tmp_path / 'wbp_custom_filter.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'wbp', '--filter', 'hamming', '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)


class TestGradientDescent:

    def test_sart_end_to_end(self, tmp_path):
        outPath = tmp_path / 'sart.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'sart', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_sirt_end_to_end(self, tmp_path):
        outPath = tmp_path / 'sirt.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'sirt', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_ossart_end_to_end(self, tmp_path):
        outPath = tmp_path / 'ossart.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'ossart', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_asdpocs_end_to_end(self, tmp_path):
        outPath = tmp_path / 'asdpocs.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'ASD-POCS', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_osasdpocs_end_to_end(self, tmp_path):
        outPath = tmp_path / 'osasdpocs.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'osasdpocs', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_pcsd_end_to_end(self, tmp_path):
        outPath = tmp_path / 'pcsd.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'pcsd', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_awpcsd_end_to_end(self, tmp_path):
        outPath = tmp_path / 'awpcsd.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'awpcsd', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_awasdpocs_end_to_end(self, tmp_path):
        outPath = tmp_path / 'awasdpocs.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'awasdpocs', '--iter', str(ITERATIONS_GRADIENT), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)


class TestKrylov:

    def test_cgls_end_to_end(self, tmp_path):
        outPath = tmp_path / 'cgls.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'cgls', '--iter', str(ITERATIONS_KRYLOV), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_lsqr_end_to_end(self, tmp_path):
        outPath = tmp_path / 'lsqr.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'lsqr', '--iter', str(ITERATIONS_KRYLOV), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_lsmr_end_to_end(self, tmp_path):
        outPath = tmp_path / 'lsmr.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'lsmr', '--iter', str(ITERATIONS_KRYLOV), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_hybridlsqr_end_to_end(self, tmp_path):
        outPath = tmp_path / 'hybridlsqr.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'hybridlsqr', '--iter', str(ITERATIONS_KRYLOV), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_abgmres_end_to_end(self, tmp_path):
        outPath = tmp_path / 'abgmres.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'abgmres', '--iter', str(ITERATIONS_KRYLOV), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_bagmres_end_to_end(self, tmp_path):
        outPath = tmp_path / 'bagmres.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'bagmres', '--iter', str(ITERATIONS_KRYLOV), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_irntvcgls_end_to_end(self, tmp_path):
        outPath = tmp_path / 'irntvcgls.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'irntvcgls', '--iter', str(ITERATIONS_KRYLOV), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)


class TestStatistical:

    def test_mlem_end_to_end(self, tmp_path):
        outPath = tmp_path / 'mlem.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'mlem', '--iter', str(ITERATIONS_STATISTICAL), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)


class TestVariational:

    def test_fista_end_to_end(self, tmp_path):
        outPath = tmp_path / 'fista.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'fista', '--iter', str(ITERATIONS_VARIATIONAL), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)

    def test_sarttv_end_to_end(self, tmp_path):
        outPath = tmp_path / 'sarttv.mrc'

        result = runCli([
            '--tiltseries', TS_PATH, '--angles', ANGLES_PATH, '--xf', XF_PATH,
            '--thickness', str(TOMO_THICKNES), '-o', str(outPath),
            '--method', 'sarttv', '--iter', str(ITERATIONS_VARIATIONAL), '--gpu', GPU_ID,
        ])

        assertReconstructionOk(result, outPath)
