from __future__ import division
import numpy as np
import numpy.matlib as matlib

class geometry:

    def __init__(self):

        self.mode = None
        self.COR = None
        self.accuracy = 0.5
    def check_geo(self,angles):

        if self.mode =='cone':
            manditory_attribs = ['nVoxel','sVoxel','dVoxel',
                                'nDetector', 'sDetector', 'dDetector',
                                'DSO','DSD']
            included_attribs_indx = [hasattr(self,attrib) for attrib in manditory_attribs]
            if not all(included_attribs_indx):
                raise AttributeError('following manditory fields '
                                     'missing from geometry:' + str([attrib for attrib in manditory_attribs
                                                                     if not hasattr(self,attrib)])
                )
            optional_fields = ['offOrigin','offDetector','rotDetector','COR',
                               'mode','accuracy']

            #image data
            if not self.nVoxel.shape == (3, ): raise AttributeError('geo.nVoxel.shape should be (3, )')
            if not self.sVoxel.shape == (3, ): raise AttributeError('geo.sVoxel.shape should be (3, )')
            if not self.dVoxel.shape == (3, ): raise AttributeError('geo.dVoxel.shape should be (3, )')
            if not sum(abs(self.dVoxel*self.nVoxel - self.sVoxel)) < 1e-6: 'nVoxel*dVoxel is not equal to sVoxel. ' \
                                                                           'Check fields.'

            #Detector Data
            if not self.nDetector.shape == (2, ): raise AttributeError('geo.nDecetor.shape should be (2, )')
            if not self.sDetector.shape == (2, ): raise AttributeError('geo.sDetector.shape should be (2, )')
            if not self.dDetector.shape == (2, ): raise AttributeError('geo.dDetector.shape should be (2, )')
            if not sum(abs(self.dDetector*self.nDetector - self.sDetector)) < 1e-6: raise AttributeError(
                'nDetector*dDetecor is not equal to sVoxel. Check fields.')

            #TODO: DSD, DSO need to be implemented once source code is updated
            for attrib in ['DSD','DSO']:
                self._check_and_repmat(attrib,angles)
            if hasattr(self,'offOrigin'):
                self._check_and_repmat('offOrigin',angles)
            if hasattr(self,'offDetector'):
                self._check_and_repmat('offDetector',angles)
            if hasattr(self, 'rotDetector'):
                self._check_and_repmat('rotDetector',angles)
            if self.COR != None:
                self._check_and_repmat('COR',angles)




        if self.mode == 'paralell':
            pass
        elif self.mode == None:
            raise AttributeError("geometry.mode needs to be specified in order for tigre to assess whether the geometry is "
                                 "correctly implemented")


    def _check_and_repmat(self, attrib, angles):
        """
        Checks whether the attribute is a single value and repeats it into an array if it is
        :param attrib: string
        :param angles: np.ndarray
        """
        old_attrib = getattr(self,attrib)

        if type(old_attrib) in [float,int,np.float32]:
            new_attrib = matlib.repmat(old_attrib,1,len(angles))[0]
            setattr(self,attrib,new_attrib)

        elif type(old_attrib) == np.ndarray:
            if len(old_attrib.shape)==1:
                new_attrib=matlib.repmat(old_attrib,len(old_attrib),len(angles))
                setattr(self,attrib,new_attrib)

            elif old_attrib.shape not in [(max(angles.shape), ),
                                          (1, max(angles.shape)),
                                          (old_attrib.shape[0], max(angles.shape))]:

                raise AttributeError(attrib + "array of shape " + str(old_attrib.shape) + " not compatible for " + attrib)

    def __str__(self):
        parameters = []
        parameters.append("TIGRE parameters")
        parameters.append("-----")
        parameters.append("Distance from source to detector = " + str(self.DSD) + "mm")
        parameters.append("Distance from source to origin = " + str(self.DSO) + "mm")

        parameters.append("-----")
        parameters.append("Detector parameters")
        parameters.append("Number of pixels = " + str(self.nDetector))
        parameters.append("Size of each pixel = " + str(self.dDetector) + "mm")
        parameters.append("Total size of the detector = " + str(self.sDetector) + "mm")

        parameters.append("-----")
        parameters.append("Image parameters")
        parameters.append("Number of voxels = " + str(self.nVoxel))
        parameters.append("Total size of the image = " + str(self.sVoxel) + "mm")
        parameters.append("Size of each voxel = " + str(self.dVoxel) + "mm")

        parameters.append("-----")
        parameters.append("Offset correction parameters")
        parameters.append("Offset of image from origin = " + str(self.offOrigin) + "mm")
        parameters.append("Offset of detector = " + str(self.offDetector) + "mm")

        parameters.append("-----")
        parameters.append("Auxillary parameters")
        parameters.append("Accuracy of forward projection = " + str(self.accuracy))

        return '\n'.join(parameters)
