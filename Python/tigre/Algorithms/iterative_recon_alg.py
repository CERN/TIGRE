import numpy as np
from data_minimization import DataMinimization
from regularisation import Regularisation

class IterativeReconAlg(Regularisation, DataMinimization):

    def __init__(self,proj,geo,angles,niter,**kwargs):

        self.proj = proj
        self.angles = angles
        self.geo = geo
        self.niter = niter

        options = dict(blocksize=20, lmbda=1, lmbda_red=0.99, OrderStrategy=None, Quameasopts=None, init=None,
                       verbose=True, noneg=True, computel2=False)
        for kw in kwargs.keys()
            if not options.has_key(kw):
                #Note: might not want this warning (typo checking).
                print("Warning: " + kw + "not recognised as default parameter for instance of IterativeReconAlg.")
        self.__dict__.update(options)
        self.__dict__.update(**kwargs)

    def get_w(self):
        pass










if __name__=='__main__':
    new_class = IterativeReconAlg('proj', 'geo', 'angles', 'niter', **{'e': 1, 'f': 2, 'g': 4})
    new_class.change_proj()
    print(new_class.__dict__)
