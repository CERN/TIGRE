import os
import unittest

import numpy as np
import tigre
from algorithm_test import AlgorithmTest

dirname = os.path.dirname(__file__)
targetdir = str(np.load(os.path.join(dirname, "targetdir.npy")))
algs = tigre.algorithms.__all__
numconfigs = 4


class TestSequence(unittest.TestCase):
    pass


def test_generator(configuration, algorithm):
    def test(self):
        if algorithm == "fbp" and (configuration in ["configuration1.npy", "configuration3.npy"]):
            unittest.skip("FBP skipped for conebeam")
        elif algorithm == "FDK" and (configuration in ["configuration2.npy", "configuration4.npy"]):
            unittest.skip("FDK skipped for parallel beam. ")
        else:
            self.assertTrue(AlgorithmTest(configuration, algorithm).unit_test_call())

    return test


if __name__ == "__main__":
    for alg in algs:
        for num in range(1, numconfigs + 1):
            test_name = "test_%d_%s" % (num, alg)
            test = test_generator("configuration" + str(num) + ".npy", alg)
            setattr(TestSequence, test_name, test)
    import xmlrunner

    unittest.main(testRunner=xmlrunner.XMLTestRunner(output=os.path.join(dirname, "test_reports")))
