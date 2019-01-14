from __future__ import division
from __future__ import print_function
import numpy as np
import numpy.matlib as matlib
import inspect


class Geometry(object):

    def __init__(self):
        self.mode = None
        self.accuray = 0.5
        self.n_proj = None
        self.angles = None
        self.filter = None
        self.rotDetector = np.array((0, 0, 0))

    def check_geo(self, angles, verbose=False):
        if angles.ndim == 1:
            self.n_proj = angles.shape[0]
            zeros_array = np.zeros((self.n_proj, 1), dtype=np.float32)
            self.angles = np.hstack((angles.reshape(self.n_proj, 1), zeros_array, zeros_array))

        elif angles.ndim == 2:
            if angles.shape[1] != 3:
                raise BufferError("Expected angles of dimensions (n, 3), got: " + str(angles.shape))
            self.n_proj = angles.shape[0]
            angles = angles.copy()
            setattr(self, 'angles', angles)
        else:
            raise BufferError("Unexpected angles shape: " + str(angles.shape))
        if self.mode == None:
            setattr(self, 'mode', 'cone')

        manditory_attribs = ['nVoxel', 'sVoxel', 'dVoxel',
                             'nDetector', 'sDetector', 'dDetector',
                             'DSO', 'DSD']
        included_attribs_indx = [hasattr(self, attrib) for attrib in manditory_attribs]
        if not all(included_attribs_indx):
            raise AttributeError('following manditory fields '
                                 'missing from geometry:' + str([attrib for attrib in manditory_attribs
                                                                 if not hasattr(self, attrib)])
                                 )
        optional_attribs = ['offOrigin', 'offDetector', 'rotDetector', 'COR',
                            'mode', 'accuracy']

        # image data
        if not self.nVoxel.shape == (3,): raise AttributeError('geo.nVoxel.shape should be (3, )')
        if not self.sVoxel.shape == (3,): raise AttributeError('geo.sVoxel.shape should be (3, )')
        if not self.dVoxel.shape == (3,): raise AttributeError('geo.dVoxel.shape should be (3, )')
        if not sum(abs(self.dVoxel * self.nVoxel - self.sVoxel)) < 1e-6: 'nVoxel*dVoxel is not equal to sVoxel. ' \
                                                                         'Check fields.'

        # Detector Data
        if not self.nDetector.shape == (2,): raise AttributeError('geo.nDecetor.shape should be (2, )')
        if not self.sDetector.shape == (2,): raise AttributeError('geo.sDetector.shape should be (2, )')
        if not self.dDetector.shape == (2,): raise AttributeError('geo.dDetector.shape should be (2, )')
        if not sum(abs(self.dDetector * self.nDetector - self.sDetector)) < 1e-6: raise AttributeError(
            'nDetector*dDetecor is not equal to sVoxel. Check fields.')

        for attrib in ['DSD', 'DSO']:
            self._check_and_repmat(attrib, angles)

        if hasattr(self, 'offOrigin'):
            self._check_and_repmat('offOrigin', angles)
        else:
            self.offOrigin = np.array([0, 0, 0])
            self._check_and_repmat('offOrigin', angles)

        if hasattr(self, 'offDetector'):
            self._check_and_repmat('offDetector', angles)
        else:
            self.offDetector = np.zeros((angles.shape[0], 2))

        if hasattr(self, 'rotDetector'):
            self._check_and_repmat('rotDetector', angles)
        else:
            self.rotDetector = np.zeros((angles.shape[0], 3))

        if hasattr(self, 'COR'):
            self._check_and_repmat('COR', angles)
        else:
            self.COR = np.zeros(angles.shape[0])

        if verbose:
            self._verbose_output()

    def _check_and_repmat(self, attrib, angles):
        """
        Checks whether the attribute is a single value and repeats it into an array if it is
        :param attrib: string
        :param angles: np.ndarray
        """
        old_attrib = getattr(self, attrib)

        if type(old_attrib) in [float, int, np.float32, np.float64]:
            new_attrib = matlib.repmat(old_attrib, 1, angles.shape[0])[0]
            setattr(self, attrib, new_attrib)

        elif type(old_attrib) == np.ndarray:
            if old_attrib.ndim == 1:
                if old_attrib.shape in [(3,), (2,), (1,)]:
                    new_attrib = matlib.repmat(old_attrib, angles.shape[0], 1)
                    setattr(self, attrib, new_attrib)
                elif old_attrib.shape == (angles.shape[0],):
                    pass
            else:
                if old_attrib.shape == (angles.shape[0], old_attrib.shape[1]):
                    pass
                else:
                    raise AttributeError(attrib + " with shape: " + str(old_attrib.shape) +
                                         " not compible with shapes: " + str([(angles.shape[0],),
                                                                              (angles.shape[0], old_attrib.shape[1]),
                                                                              (3,), (2,), (1,)]))

        else:
            TypeError(
                "Data type not understood for: geo." + attrib + " with type = " + str(type(getattr(self, attrib))))

    def _verbose_output(self):
        for obj in inspect.getmembers(self):
            if obj[0][0] == '_':
                pass
            elif obj[0] == 'check_geo':
                pass
            elif type(obj[1]) == np.ndarray:
                print(self.mode + ': ' + str((obj[0], obj[1].shape)))
            else:
                print(self.mode + ': ' + str(obj))

    def convert_contig_mode(self):
        dim_attribs = ['nVoxel', 'sVoxel', 'dVoxel',
                       'nDetector', 'sDetector', 'dDetector']
        for attrib in dim_attribs:
            setattr(self, attrib, getattr(self, attrib)[::-1].copy())

    def issame(self, other_geo_dict):
        geo_comp_list = []
        if set(other_geo_dict) == set(self.__dict__):
            for key in self.__dict__.keys():
                if type(self.__dict__[key]) == np.ndarray:
                    geo_comp_list.append(all(self.__dict__[key] == other_geo_dict[key]))
                else:
                    geo_comp_list.append(self.__dict__[key] == other_geo_dict[key])
            return all(geo_comp_list)
        else:
            return False

    def __str__(self):
        parameters = []
        parameters.append("TIGRE parameters")
        parameters.append("-----")
        parameters.append("Distance from source to detector = " + str(self.DSD) + " mm")
        parameters.append("Distance from source to origin = " + str(self.DSO) + " mm")

        parameters.append("-----")
        parameters.append("Detector parameters")
        parameters.append("Number of pixels = " + str(self.nDetector))
        parameters.append("Size of each pixel = " + str(self.dDetector) + " mm")
        parameters.append("Total size of the detector = " + str(self.sDetector) + " mm")

        parameters.append("-----")
        parameters.append("Image parameters")
        parameters.append("Number of voxels = " + str(self.nVoxel))
        parameters.append("Total size of the image = " + str(self.sVoxel) + " mm")
        parameters.append("Size of each voxel = " + str(self.dVoxel) + " mm")

        parameters.append("-----")
        parameters.append("Offset correction parameters")
        parameters.append("Offset of image from origin = " + str(self.offOrigin) + " mm")
        parameters.append("Offset of detector = " + str(self.offDetector) + " mm")

        parameters.append("-----")
        parameters.append("Auxillary parameters")
        parameters.append("Accuracy of forward projection = " + str(self.accuracy))

        return '\n'.join(parameters)


class ParallelGeo(Geometry):
    def __init__(self,nVoxel):
        if nVoxel is None:
            raise ValueError('nVoxel needs to be given for initialisation of parallel beam')
        Geometry.__init__(self)
        self.nVoxel = nVoxel
        self.dVoxel = np.array([1,1,1],dtype=np.float32)
        self.dDetector = np.array([1,1],dtype=np.float32)
        self.DSO = self.nVoxel[0]
        self.DSD = self.nVoxel[0]*2
        self.nDetector = self.nVoxel[:2]
        self.sVoxel = self.nVoxel
        self.sDetector = self.nVoxel[:2]
        self.accuracy = 0.5


def geometry(mode='cone',nVoxel=None):
    if mode in [None, 'cone']:
        return Geometry()
    else:
        return ParallelGeo(nVoxel)
