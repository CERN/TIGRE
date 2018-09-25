from __future__ import division
from matplotlib import pyplot as plt
import numpy as np


class plotImg:
    # NOTE: type help(plotImg) after importing in order to get a readable manual.
    ('\n'
     'plotImg(cube, dim) \n'
     '    plots figure \n'
     'default: progressive in slices following\n'
     '    axis (dim)\n'
     'Parameters \n'
     '---------- \n'
     'cube : Any 3D numpy array \n'
     '\n'
     'dim : ("X","Y","Z","x","y","z"), optional \n'
     '       default is "Z"\n'
     '       NOTE: string arguments!'
     '\n'
     'slice: int, optional\n'
     '     returns page of matrix according to index\n'
     'Examples:\n'
     '---------\n'
     'a=np.ones([3,3,3])\n'
     'plotImg(a)\n'
     '>>>returns plot along dim Z\n'
     'plotImg(a,dim="X")\n'
     '>>>returns plot along dim X\n')

    def __init__(self, cube, dim=None, slice=None):
        self.cube = cube
        self.dim = dim
        self.slice = slice
        self.dimint = None  # keeps track of what dim
        self.dimlist = ['X', 'Y', 'Z', 'x', 'y', 'z', None]  # accepted parameters for dim
        if self.slice is None:
            self.run()
        if self.slice is not None:
            self.slicer()

    def run(self):
        if self.dim not in self.dimlist:
            raise NameError('check inputs for dim, should be string.')

        if self.dim in [None, 'Z', 'z']:
            self.dimint = 2
            self.dimlist = ['X', 'Y', 'Z']
            self.run_plot()
        if self.dim in ['X', 'x']:
            self.dimint = 0
            self.dimlist = ['Y', 'Z', 'X']
            self.run_plot()
        if self.dim in ['Y', 'y']:
            self.dimint = 1
            self.dimlist = ['X', 'Z', 'Y']
            self.run_plot()

    def run_plot(self):
        plt.ion()
        min_val = np.amin(self.cube)
        max_val = np.amax(self.cube)
        dim = self.cube.shape

        for i in range(dim[self.dimint]):
            plt.clf()
            if self.dimint == 2:
                plt.imshow(np.squeeze(self.cube[:, :, i]), cmap=plt.cm.gray, origin='lower', vmin=min_val, vmax=max_val)
            if self.dimint == 1:
                plt.imshow(np.squeeze(self.cube[:, i].transpose()), cmap=plt.cm.gray, origin='lower', vmin=min_val,
                           vmax=max_val)
            if self.dimint == 0:
                plt.imshow(np.squeeze(self.cube[i].transpose()), cmap=plt.cm.gray, origin='lower', vmin=min_val,
                           vmax=max_val)
            plt.gca().get_xaxis().set_ticks([])
            plt.gca().get_yaxis().set_ticks([])
            plt.xlabel(self.dimlist[0])
            plt.ylabel(self.dimlist[1])
            plt.title(self.dimlist[2] + ':' + str(i))

            plt.colorbar()

            plt.pause(0.01)

    def slicer(self):
        # NOTE: Transpose not quite right (still returning image at wrong angle, No labels for axes)
        min_val = np.amin(self.cube)
        max_val = np.amax(self.cube)
        if self.dim in [None, 'Z', 'z']:
            plt.imshow(np.squeeze(self.cube[:,:, self.slice]), cmap=plt.cm.gray, origin='lower', vmin=min_val, vmax=max_val)
        if self.dim in ['Y', 'y']:
            plt.imshow(np.squeeze(self.cube[:, self.slice]), cmap=plt.cm.gray,origin='lower', vmin=min_val, vmax=max_val)
        if self.dim in ['X', 'x']:
            plt.imshow(np.squeeze(self.cube[self.slice]).transpose(), cmap=plt.cm.gray,origin='lower', vmin=min_val, vmax=max_val)
        plt.show()
