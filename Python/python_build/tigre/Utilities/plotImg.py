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
     'dim : int, describes axis of progression.'
     '\n'
     'slice: int, optional\n'
     '     returns page of matrix according to index\n'
     'slider: bool, optional\n'
     '     returns interactive plot'
     'Examples:\n'
     '---------\n'
     'a=np.ones([3,3,3])\n'
     'plotImg(a)\n'
     '>>>returns plot along dim 0\n'
     'plotImg(a,dim=0)\n'
     '>>>returns plot along dim X\n')

    def __init__(self, cube, dim=None, slice=None,slider=False,**kwargs):
        self.cube = cube
        self.dim = dim
        self.slice = slice
        self.slider=slider
        self.dimint = dim
        self.kwargs=kwargs
        if self.slider:
            self.cube_show_slider()
            return 
        if self.slice is None:
            self.run_plot()
        if self.slice is not None:
            self.slicer()
    def run_plot(self):
        if self.dimint is None:
            self.dimint = 0

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
            plt.title('slice' + ':' + str(i))

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
        kwargs.update(self.kwargs)
        cube=self.cube
        if self.dimint is None:
            axis=0
        else:
            axis=self.dimint


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

