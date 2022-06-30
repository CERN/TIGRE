from __future__ import division
from __future__ import print_function

import inspect

import numpy as np
import numpy.matlib as matlib
import tigre


class Geometry(object):
    def __init__(self):
        self.mode = None
        self.n_proj = None
        self.angles = None
        self.filter = None
        self.nVoxel = np.zeros(3)
        self.dVoxel = np.zeros(3)
        self.sVoxel = np.zeros(3)
        self.nDetector = np.zeros(2)
        self.dDetector = np.zeros(2)
        self.sDetector = np.zeros(2)


    def check_geo(self, angles, verbose=False):
        if angles.ndim == 1:
            self.n_proj = angles.shape[0]
            zeros_array = np.zeros((self.n_proj, 1), dtype=np.float32)
            self.angles = np.hstack((angles.reshape(self.n_proj, 1), zeros_array, zeros_array))

        elif angles.ndim == 2:
            if angles.shape[1] != 3:
                raise BufferError("Expected angles of dimensions (n, 3), got: " + str(angles.shape))
            self.n_proj = angles.shape[0]
            self.angles = angles.copy()
        else:
            raise BufferError("Unexpected angles shape: " + str(angles.shape))
        if self.mode is None:
            self.mode = "cone"

        mandatory_attribs = [
            "nVoxel",
            "sVoxel",
            "dVoxel",
            "nDetector",
            "sDetector",
            "dDetector",
            "DSO",
            "DSD",
        ]
        included_attribs_indx = [hasattr(self, attrib) for attrib in mandatory_attribs]
        if not all(included_attribs_indx):
            raise AttributeError(
                "following manditory fields "
                "missing from geometry:"
                + str([attrib for attrib in mandatory_attribs if not hasattr(self, attrib)])
            )
        optional_attribs = ["offOrigin", "offDetector", "rotDetector", "COR", "mode", "accuracy"]

        # image data
        if not self.nVoxel.shape == (3,):
            raise AttributeError("geo.nVoxel.shape should be (3, )")
        if not self.sVoxel.shape == (3,):
            raise AttributeError("geo.sVoxel.shape should be (3, )")
        if not self.dVoxel.shape == (3,):
            raise AttributeError("geo.dVoxel.shape should be (3, )")
        if not sum(abs(self.dVoxel * self.nVoxel - self.sVoxel)) < 1e-6:
            "nVoxel*dVoxel is not equal to sVoxel. " "Check fields."

        # Detector Data
        if not self.nDetector.shape == (2,):
            raise AttributeError("geo.nDecetor.shape should be (2, )")
        if not self.sDetector.shape == (2,):
            raise AttributeError("geo.sDetector.shape should be (2, )")
        if not self.dDetector.shape == (2,):
            raise AttributeError("geo.dDetector.shape should be (2, )")
        if not sum(abs(self.dDetector * self.nDetector - self.sDetector)) < 1e-6:
            raise AttributeError("nDetector*dDetecor is not equal to sDetector. Check fields.")

        for attrib in ["DSD", "DSO"]:
            self.__check_and_repmat__(attrib, angles)

        if hasattr(self, "offOrigin"):
            self.__check_and_repmat__("offOrigin", angles)
        else:
            self.offOrigin = np.array([0, 0, 0])
            self.__check_and_repmat__("offOrigin", angles)

        if hasattr(self, "offDetector"):
            self.__check_and_repmat__("offDetector", angles)
        else:
            self.offDetector = np.array([0, 0])
            self.__check_and_repmat__("offDetector", angles)

        if hasattr(self, "rotDetector"):
            self.__check_and_repmat__("rotDetector", angles)
        else:
            self.rotDetector = np.array([0, 0, 0])
            self.__check_and_repmat__("rotDetector", angles)

        if hasattr(self, "COR"):
            self.__check_and_repmat__("COR", angles)
        else:
            self.COR = np.zeros(angles.shape[0])

        if not hasattr(self, "accuracy"):
            self.accuracy = 0.5
        # IMPORTANT: cast all numbers to float32
        if verbose:
            self._verbose_output()

    def checknans(self):
        for attrib in self.__dict__:
            if str(getattr(self, attrib)) == "nan":
                raise ValueError("nan found for Geometry abbtribute:" + attrib)
            elif type(getattr(self, attrib)) == np.ndarray:
                if np.isnan(getattr(self, attrib)).all():
                    raise ValueError("Nan found in Geometry abbtribute:" + attrib)

    def cast_to_single(self):
        """
        Casts all number values in current instance to single prevision floating point types.

        :return: None
        """
        for attrib in self.__dict__:
            if getattr(self, attrib) is not None:
                try:
                    setattr(self, attrib, np.float32(getattr(self, attrib)))
                except ValueError:
                    pass

    def __check_and_repmat__(self, attrib, angles):
        """
        Checks whether the attribute is a single value and repeats it into an array if it is

        :rtype: None
        :param attrib: string
        :param angles: np.ndarray
        """
        old_attrib = getattr(self, attrib)

        if type(old_attrib) in [float, int, np.float32, np.float64, np.int32]:
            new_attrib = matlib.repmat(old_attrib, 1, angles.shape[0])[0]
            setattr(self, attrib, new_attrib)

        elif type(old_attrib) == np.ndarray:
            if old_attrib.ndim == 1:
                if old_attrib.shape in [(3,), (2,)] and attrib not in ["DSD", "DSO", "COR"]:
                    new_attrib = matlib.repmat(old_attrib, angles.shape[0], 1)
                    setattr(self, attrib, new_attrib)
                elif old_attrib.shape in [(1,)]:
                    new_attrib = matlib.repmat(old_attrib, angles.shape[0], 1)
                    setattr(self, attrib, new_attrib)
                elif old_attrib.shape == (angles.shape[0],):
                    pass
            else:
                if old_attrib.shape == (angles.shape[0], old_attrib.shape[1]):
                    pass
                else:
                    raise AttributeError(
                        attrib
                        + " with shape: "
                        + str(old_attrib.shape)
                        + " not compatible with shapes: "
                        + str(
                            [
                                (angles.shape[0],),
                                (angles.shape[0], old_attrib.shape[1]),
                                (3,),
                                (2,),
                                (1,),
                            ]
                        )
                    )

        else:
            raise TypeError(
                "Data type not understood for: geo."
                + attrib
                + " with type = "
                + str(type(getattr(self, attrib)))
            )

    def _verbose_output(self):
        for obj in inspect.getmembers(self):
            if obj[0][0] == "_":
                pass
            elif obj[0] == "check_geo":
                pass
            elif type(obj[1]) == np.ndarray:
                print(self.mode + ": " + str((obj[0], obj[1].shape)))
            else:
                print(self.mode + ": " + str(obj))

    def __str__(self):
        parameters = []
        parameters.append("TIGRE parameters")
        parameters.append("-----")
        parameters.append("Geometry parameters")
        parameters.append("Distance from source to detector (DSD) = " + str(self.DSD) + " mm")
        parameters.append("Distance from source to origin (DSO)= " + str(self.DSO) + " mm")
        parameters.append("-----")
        parameters.append("Detector parameters")
        parameters.append("Number of pixels (nDetector) = " + str(self.nDetector))
        parameters.append("Size of each pixel (dDetector) = " + str(self.dDetector) + " mm")
        parameters.append("Total size of the detector (sDetector) = " + str(self.sDetector) + " mm")
        parameters.append("-----")
        parameters.append("Image parameters")
        parameters.append("Number of voxels (nVoxel) = " + str(self.nVoxel))
        parameters.append("Total size of the image (sVoxel) = " + str(self.sVoxel) + " mm")
        parameters.append("Size of each voxel (dVoxel) = " + str(self.dVoxel) + " mm")

        parameters.append("-----")
        if hasattr(self, "offOrigin") and hasattr(self, "offDetector"):
            parameters.append("Offset correction parameters")
            if hasattr(self, "offOrigin"):
                parameters.append(
                    "Offset of image from origin (offOrigin) = " + str(self.offOrigin) + " mm"
                )
            if hasattr(self, "offDetector"):
                parameters.append(
                    "Offset of detector (offDetector) = " + str(self.offDetector) + " mm"
                )

        parameters.append("-----")
        parameters.append("Auxillary parameters")
        parameters.append(
            "Samples per pixel of forward projection (accuracy) = " + str(self.accuracy)
        )

        if hasattr(self, "rotDetector"):
            parameters.append("-----")
            parameters.append(
                "Rotation of the Detector (rotDetector) = " + str(self.rotDetector) + " rad"
            )

        if hasattr(self, "COR"):
            parameters.append("-----")
            parameters.append("Centre of rotation correction (COR) = " + str(self.COR) + " mm")

        return "\n".join(parameters)

    def __cmp__(self, other):
        resultofnumpiesanallyretentiveattemptatbeingphilosophical = []
        for attrib in self.__dict__:
            result = getattr(self, attrib) == getattr(other, attrib)
            try:
                resultofnumpiesanallyretentiveattemptatbeingphilosophical.append(result.all())
            except Exception:
                try:
                    resultofnumpiesanallyretentiveattemptatbeingphilosophical.extend(result)
                except Exception:
                    resultofnumpiesanallyretentiveattemptatbeingphilosophical.append(result)

        # why is this boolean reversed when returned?
        # because for some reason its reversed when i return it from this function. Who knows.

        return not all(resultofnumpiesanallyretentiveattemptatbeingphilosophical)


class ParallelGeo(Geometry):
    def __init__(self, nVoxel):  # noqa: N803
        if nVoxel is None:
            raise ValueError("nVoxel needs to be given for initialisation of parallel beam")
        Geometry.__init__(self)
        self.mode = "parallel"

        self.nVoxel = nVoxel
        self.dVoxel = np.array([1, 1, 1])
        self.sVoxel = self.nVoxel

        self.DSO = np.float32(self.nVoxel[0])
        self.DSD = np.float32(self.nVoxel[0] * 2)

        self.dDetector = np.array([1, 1])
        self.nDetector = np.zeros((2,), dtype=np.int64)
        self.nDetector[0] = self.nVoxel[0]
        self.nDetector[1] = max(self.nVoxel[1], self.nVoxel[2])
        # self.nDetector = self.nVoxel[:2]
        self.sDetector = self.nDetector

        self.accuracy = 0.5

        self.offOrigin = np.array([0, 0, 0])
        self.offDetector = np.array([0, 0])

        self.rotDetector = np.array([0, 0, 0])


def geometry(mode="cone", nVoxel=None, default=False, high_resolution=True):
    """
    Constructor for geometry used in reconstruction of images in TIGRE

    Parameters
    ----------
    :param mode: (str)
        'cone' or 'parallel'
    :param nVoxel: (np.ndarray)
        number of voxels the reconstruction is composed of
    :param default: (bool)
        calculates other parameters in geometry. is by default true for
        parallel geometry
    :param high_resolution: (bool)
        preset values for geometry in mode=cone. WARNING: for smaller
        tests it is better to use this rather than setting nVoxel
        manually.
        if true: nVoxel = np.array([512,512,512])
        if false: nVoxel = np.array([64,64,64])
    :return: (tigre.geometry.Geometry)

    Usage
    -----
    >>> import tigre
    >>> #Cone beam with no preset parameters
    >>> geo_cone = tigre.geometry(mode='cone')
    >>> # Cone beam default, low quality
    >>> geo_cone_default = tigre.geometry(mode='cone',high_resolution=False)
    >>> # Cone beam with specific nVoxel requirements
    >>> geo_cone__default2 = tigre.geometry(nVoxel=np.array([64,64,64]),
    >>>                                     mode='cone'
    >>>                                     default=True)
    >>> # Parallel beam
    >>> geo_par = tigre.geometry(mode='parallel',
    >>>                          nVoxel=np.array([64,64,64]))
    """
    if mode == "cone":
        if default:
            return tigre.geometry_default(high_resolution, nVoxel)
        else:
            return Geometry()
    if mode == "parallel":
        return ParallelGeo(nVoxel)
    else:
        raise ValueError("mode: " + mode + " not recognised.")
