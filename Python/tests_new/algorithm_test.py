
import os
import numpy as np
import sys
import tigre.algorithms as algs

class AlgorithmTest(object):

    def __init__(self, configuration,algorithm,**kwargs):
        """

        :param configuration: (str)
            which predefined Configuration to use
        :param algorithm: (str)
            which algorithm to test

        """
        self.dirname = os.path.dirname(__file__)
        self.targetdir = str(np.load(os.path.join(self.dirname,'targetdir.npy')))
        configdict = np.load(os.path.join(self.dirname,configuration)).item()
        for key in configdict:
            setattr(self,key,configdict[key])
        self.algorithm = algorithm
        self.testpassed = False
        self.confignumber = os.path.splitext(configuration)[0]
    def test(self):
        if self.algorithm == 'fbp' and self.geo.mode != 'parallel':
            print('WARNING: fbp was implemented in a cone beam regime. \n')
            print('Test ignored.')
            raise SystemExit()
        else:
            print(self.geo.mode)
    def save_output(self):
        resultfilename = self.confignumber + '.npy'
        try:
            resultsdata = np.load(os.path.join(self.targetdir, resultfilename)).item()

        except Exception:
            resultsdata = dict()
        resultsdata.update({self.algorithm : self.testpassed})
        np.save(os.path.join(self.targetdir,resultfilename),resultsdata)



if __name__ == '__main__':
    configuration = sys.argv[1]
    algorithm = sys.argv[2]

    test = AlgorithmTest(configuration,algorithm)
    try:
        test.test()
    except Exception as e:
        print(e +'\n')
        print(configuration)
        print(algorithm)

    test.save_output()