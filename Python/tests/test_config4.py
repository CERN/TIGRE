import unittest
import sys
import os
import numpy as np
dirname = os.path.dirname(__file__)
targetdir = str(np.load(os.path.join(dirname, 'targetdir.npy')))


class config4TestCase(unittest.TestCase):
    def setUp(self):

        self.results = np.load(
            os.path.join(
                targetdir,
                'configuration4.npy'),
            allow_pickle=True).item()

    def testFDK(self):
        self.assertTrue(self.results['FDK'])

    def testfbp(self):
        self.assertTrue(self.results['fbp'])

    def testsirt(self):
        self.assertTrue(self.results['sirt'])

    def testossart(self):
        self.assertTrue(self.results['ossart'])

    def testcgls(self):
        self.assertTrue(self.results['cgls'])

    def testawasd_pocs(self):
        self.assertTrue(self.results['awasd_pocs'])

    def testasd_pocs(self):
        self.assertTrue(self.results['asd_pocs'])

    def testfista(self):
        self.assertTrue(self.results['fista'])


if __name__ == '__main__':
    import xmlrunner
    unittest.main(
        testRunner=xmlrunner.XMLTestRunner(
            output=os.path.join(
                dirname,
                'test_reports')))
