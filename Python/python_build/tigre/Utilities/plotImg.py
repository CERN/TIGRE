from __future__ import division
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
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

    def __init__(self, cube, dim=None, slice=None,slider=False):
        self.cube = cube
        self.dim = dim
        self.slice = slice
        self.slider=slider
        self.dimint = None  # keeps track of what dim
        self.dimlist = ['X', 'Y', 'Z', 'x', 'y', 'z', None]  # accepted parameters for dim
        if self.slider:
            self.cube_show_slider()
            return 
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

    def cube_show_slider(self, **kwargs):
        """
        Display a 3d ndarray with a slider to move along the third dimension.

        Extra keyword arguments are passed to imshow
        """
        cube=self.cube
        axis=0

        # check dim
        if not cube.ndim == 3:
            raise ValueError("cube should be an ndarray with ndim == 3")

        # generate figure
        fig = plt.figure()
        ax = plt.subplot(111)
        fig.subplots_adjust(left=0.25, bottom=0.25)

        # select first image
        s = [slice(0, 1) if i == axis else slice(None) for i in xrange(3)]
        im = cube[s].squeeze()

        # display image
        min_val = np.amin(cube)
        max_val = np.amax(cube)
        l = ax.imshow(im, cmap=plt.cm.gray, origin='lower', vmin=min_val, vmax=max_val, **kwargs)

        # define slider
        axcolor = 'lightgoldenrodyellow'
        ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

        slider = Slider(ax, 'Axis %i index' % axis, 0, cube.shape[axis] - 1,
                        valinit=0, valfmt='%i')

        def update(val):
            ind = int(slider.val)
            s = [slice(ind, ind + 1) if i == axis else slice(None)
                 for i in xrange(3)]
            im = cube[s].squeeze()
            l.set_data(im, **kwargs)
            fig.canvas.draw()

        slider.on_changed(update)

        plt.show()

