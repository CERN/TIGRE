
import os
import numpy as np
import sys
import tigre
import tigre.algorithms as algs
from tigre.demos.Test_data.data_loader import load_head_phantom
from tigre.utilities.Measure_Quality import Measure_Quality
import traceback
from matplotlib import pyplot as plt
import time


class AlgorithmTest(object):

    def __init__(self, configuration, algorithm, **kwargs):
        """

        :param configuration: (str)
            which predefined Configuration to use
        :param algorithm: (str)
            which algorithm to test

        """
        self.dirname = os.path.dirname(__file__)
        self.targetdir = str(
            np.load(
                os.path.join(
                    self.dirname,
                    'targetdir.npy'),
                allow_pickle=True))
        configdict = np.load(
            os.path.join(
                self.dirname,
                configuration),
            allow_pickle=True).item()
        for key in configdict:
            """contains: [nproj,geo,angles,niter,kwargs]"""
            setattr(self, key, configdict[key])
        self.algorithm = algorithm
        self.testpassed = False
        self.algorithm_finished = False
        self.rmse = 1.
        self.confignumber = os.path.splitext(configuration)[0]
        self.output = None
        self.timestarted = time.asctime()
        self.timeended = time.asctime()

    def test(self):
        if self.algorithm == 'fbp' and self.geo.mode != 'parallel':
            print('WARNING: fbp was implemented in cone beam.')
            print('Test ignored.\n')
            raise SystemExit()

        head = load_head_phantom(self.geo.nVoxel)
        proj = tigre.Ax(head, self.geo, self.angles)
        if self.algorithm in ['FDK', 'fbp']:
            self.output = getattr(
                tigre.algorithms,
                self.algorithm)(
                proj,
                self.geo,
                self.angles)
            self.rmse = Measure_Quality(self.output, head, ['nRMSE'])
            self.algorithm_finished = True
            return
        self.timestarted = time.asctime()
        self.output = getattr(
            tigre.algorithms,
            self.algorithm)(
            proj,
            self.geo,
            self.angles,
            self.niter,
            **self.kwargs)
        self.timeended = time.asctime()
        self.algorithm_finished = True
        self.rmse = Measure_Quality(self.output, head, ['nRMSE'])

    def unit_test_call(self):
        self.test()
        return all(self.algorithm_finished, self.rmse)

    def compound_results(self, verbose=True):
        if self.algorithm_finished and self.rmse < 0.2:
            if verbose:
                print('------------------------------------------------\n')
                print('TEST PASSED')
                print('------------------------------------------------\n')
            self.testpassed = True
        else:
            print('------------------------------------------------\n')
            print('TEST FAILED')
            print('Algorithm: ' + self.algorithm)
            print('Algorithm ran: ' + str(self.algorithm_finished))
            print('configuration number: ' + str(self.confignumber))
            print('RMSE:' + str(self.rmse))
            print('------------------------------------------------\n')

    def save_output(self):
        resultfilename = self.confignumber + '.npy'
        try:
            resultsdata = np.load(
                os.path.join(
                    self.targetdir,
                    resultfilename),
                allow_pickle=True).item()

        except Exception:
            resultsdata = dict()
        resultsdata.update({self.algorithm: self.testpassed})
        np.save(os.path.join(self.targetdir, resultfilename), resultsdata)
        if not self.testpassed:
            self.write_to_log()

    def save_fig(self):
        res = self.output
        geo = self.geo
        plt.figure()
        plt.subplot(3, 1, 1)
        plt.imshow(res[geo.nVoxel[0] / 2])
        plt.title('results for ' + self.algorithm)
        plt.ylabel('dim 0')

        plt.subplot(3, 1, 2)
        plt.imshow(res[:, geo.nVoxel[1] / 2])
        plt.ylabel('dim 1')

        plt.subplot(3, 1, 3)
        plt.imshow(res[:, :, geo.nVoxel[2] / 2])
        plt.ylabel('dim 2')
        plt.savefig(
            os.path.join(
                self.targetdir,
                self.algorithm +
                self.confignumber))

    def write_to_log(self):
        configlogfile = self.confignumber + '.log'
        logflist = []
        if configlogfile not in os.listdir(self.targetdir):
            logflist.append('GEOMETRY used for instance of testandlog: \n')
            for item in self.geo.__dict__:
                logflist.append(
                    item + ': ' + str(getattr(self.geo, item)) + '\n')
            logflist.append('nproj: ' +
                            str(self.angles.shape[0]) +
                            ' niter: ' +
                            str(self.niter) +
                            '\n')
            logflist.append(
                '------------------------------------------------\n')

        else:
            logflist.extend(
                open(
                    os.path.join(
                        self.targetdir,
                        configlogfile),
                    'r').readlines())
        logflist.append(str(self.algorithm).upper() +
                        ' ' + str(self.timestarted) + '\n')
        logflist.append('RMSE: ' + str(self.rmse) + '\n')
        logflist.append('Algorithm ran: ' +
                        str(self.algorithm_finished) + '\n')
        if self.algorithm_finished:
            logflist.append('ENDED: ' + str(self.timeended) + '\n')
        logflist.append('------------------------------------------------\n')
        logf = open(os.path.join(self.targetdir, configlogfile), 'w')
        logf.write(''.join(logflist))
        logf.close()


if __name__ == '__main__':
    configuration = sys.argv[1]
    algorithm = sys.argv[2]

    test = AlgorithmTest(configuration, algorithm)
    try:
        test.test()
    except Exception as e:
        formatedlines = traceback.format_exc()
        print(formatedlines)
    test.compound_results()
    test.save_output()
