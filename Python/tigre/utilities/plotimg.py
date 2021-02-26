from __future__ import division

import matplotlib.animation as animation
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class plotImg:
    """
    plotImg(cube, dim) 
        plots figure 
    default: progressive in slices following
        axis (dim)

    Parameters
    ----------
    cube : Any 3D numpy array

    dim : ("X","Y","Z","x","y","z"), optional
        default is "X"
        NOTE: string arguments!

    slice: int, optional
        returns page of matrix according to index

    step: int, optional
        Sets the step size between slice and slice.
        Step is 1 by default.

    savegif: string, optional
            Saves the image as .gif with the file name

    Examples:
    ---------
    a=np.ones([3,3,3])
    plotImg(a)
    >>>returns plot along dim Z

    plotImg(a,dim="X")
    >>>returns plot along dim X
    """

    def __init__(self, cube, dim=None, slice=None, step=1, savegif=None):
        self.cube = cube
        self.dim = dim
        self.slice = slice
        self.dimint = None  # keeps track of what dim
        self.dimlist = ['X', 'Y', 'Z', 'x', 'y', 'z', None]  # accepted parameters for dim
        self.step = step
        self.savegif = savegif
        if self.step is None or self.step==0:
            self.step=1
        if self.savegif=='':
            self.savegif == None
        if self.slice is None:
            self.run()
        if self.slice is not None:
            self.slicer()

    def run(self):
        if self.dim not in self.dimlist:
            raise NameError("check inputs for dim, should be string.")

        if self.dim in [None, "X", "x"]:
            self.dimint = 2
            self.dimlist = ["->Y", "->Z", "X"]
            self.run_plot()
        if self.dim in ["Y", "y"]:
            self.dimint = 1
            self.dimlist = ["->X", "->Z", "Y"]
            self.run_plot()
        if self.dim in ["Z", "z"]:
            self.dimint = 0
            self.dimlist = ["->X", "->Y", "Z"]
            self.run_plot()

    def update_frame(self, it, fig, min_val, max_val):
        i = range(0, self.cube.shape[self.dimint])[:: self.step][it]
        fig.clf()
        axis = fig.add_subplot(1, 1, 1)
        if self.dimint == 2:
            mappable = axis.imshow(
                np.squeeze(self.cube[:, :, i]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dimint == 1:
            mappable = axis.imshow(
                np.squeeze(self.cube[:, i]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dimint == 0:
            mappable = axis.imshow(
                np.squeeze(self.cube[i]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        axis.get_xaxis().set_ticks([])
        axis.get_yaxis().set_ticks([])
        axis.set_xlabel(self.dimlist[0])
        axis.set_ylabel(self.dimlist[1])
        axis.set_title(self.dimlist[2] + ":" + str(i))
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(mappable, cax=cax)
        # plt.pause(0.01)

    def run_plot(self):
        min_val = np.amin(self.cube)
        max_val = np.amax(self.cube)
        dim = self.cube.shape

        fig = plt.figure()
        ani = animation.FuncAnimation(
            fig,
            self.update_frame,
            fargs=(fig, min_val, max_val),
            interval=100,
            repeat_delay=1000,
            frames=len(range(0, dim[self.dimint])[:: self.step]),
        )
        if self.savegif is not None:
            ani.save(self.savegif, writer="pillow")
            plt.show()
        else:
            plt.show()

    def slicer(self):
        min_val = np.amin(self.cube)
        max_val = np.amax(self.cube)
        if self.dim in [None, "X", "x"]:
            plt.xlabel("Y")
            plt.ylabel("Z")
            plt.imshow(
                np.squeeze(self.cube[:, :, self.slice]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dim in ["Y", "y"]:
            plt.xlabel("X")
            plt.ylabel("Z")
            plt.imshow(
                np.squeeze(self.cube[:, self.slice]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dim in ["Z", "z"]:
            plt.xlabel("X")
            plt.ylabel("Y")
            plt.imshow(
                np.squeeze(self.cube[self.slice]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        plt.show()


plotimg = plotImg
