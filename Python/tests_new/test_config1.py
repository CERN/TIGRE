import unittest
import sys
import os
import numpy as np
dirname = os.path.dirname(__file__)
targetdir = str(np.load(os.path.join(dirname,'targetdir.npy')))

class config1TestCase(unittest.TestCase):
    def setUp(self):

        self.results = np.load(os.path.join(targetdir,'configuration1.npy')).item()
    def testFDK(self):
        self.assertTrue(self.results['FDK'])






if __name__ == '__main__':
    import xmlrunner
    unittest.main(testRunner=xmlrunner.XMLTestRunner(output=os.path.join(dirname,'test_reports')))