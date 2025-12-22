import StringIO
import os
import sys
import unittest

import numpy as np
import tigre
import tigre.algorithms as algs
from tigre.demos.Test_data import data_loader

dirname = os.path.dirname(__file__)


class TestStdout(unittest.TestCase):
    pass


def test_generator(algorithm, proj, geo, angles, niter):
    def test(self):
        capturedOutput = StringIO.StringIO()
        sys.stdout = capturedOutput
        getattr(algs, algorithm)(proj, geo, angles, niter=niter, verbose=False)
        self.assertIs(capturedOutput.getvalue(), "")
        sys.stdout = sys.__stdout__

    return test


if __name__ == "__main__":
    geo = tigre.geometry(mode="cone", default=True, high_quality=False)
    print(geo)

    true_img = data_loader.load_head_phantom(geo.nVoxel)
    angles = np.linspace(0, 2 * np.pi, 100)
    niter = 5
    proj = tigre.Ax(true_img, geo, angles)
    for alg in algs.__all__:
        if alg != "fbp":
            test_name = "test_print_%s" % (alg)
            test = test_generator(alg, proj, geo, angles, niter)
            setattr(TestStdout, test_name, test)

    unittest.main()
