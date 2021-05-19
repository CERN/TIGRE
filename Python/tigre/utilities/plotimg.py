from __future__ import division

import matplotlib
import matplotlib.animation as animation
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class plotImg:  # noqa: N801
    """
    plotImg(cube, dim)
        plots figure
    default: progressive in slices following
        axis (dim)
    """

    def __init__(
        self,
        cube,
        dim=None,
        slice=None,
        step=1,
        savegif=None,
        colormap="gray",
        clims=None,
        show_plot=None,
    ):
        self.cube = cube
        self.dim = dim
        self.slice = slice
        self.dimint = None  # keeps track of what dim
        self.dimlist = ["X", "Y", "Z", "x", "y", "z", None]  # accepted parameters for dim
        self.step = step
        self.savegif = savegif
        self.colormap = colormap
        if clims is None:
            self.min_val = np.amin(self.cube)
            self.max_val = np.amax(self.cube)
        else:
            self.min_val = clims[0]
            self.max_val = clims[1]
        if show_plot is None:
            # https://matplotlib.org/stable/tutorials/introductory/usage.html#backends
            backend = matplotlib.get_backend()
            if backend in [
                "GTK3Agg",
                "GTK3Cairo",
                "MacOSX",
                "nbAgg",
                "Qt4Agg",
                "Qt4Cairo",
                "Qt5Agg",
                "Qt5Cairo",
                "TkAgg",
                "TkCairo",
                "WebAgg",
                "WX",
                "WXAgg",
                "WXCairo",
                "module://ipykernel.pylab.backend_inline",
            ]:
                self.show_plot = True
            elif backend in ["agg", "cairo", "pdf", "pgf", "ps", "svg", "template"]:
                self.show_plot = False
            else:
                self.show_plot = True
        else:
            self.show_plot = show_plot
        if self.step is None or self.step == 0:
            self.step = 1
        if self.savegif == "":
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
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dimint == 1:
            mappable = axis.imshow(
                np.squeeze(self.cube[:, i]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dimint == 0:
            mappable = axis.imshow(
                np.squeeze(self.cube[i]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
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

        dim = self.cube.shape

        fig = plt.figure()
        ani = animation.FuncAnimation(
            fig,
            self.update_frame,
            fargs=(fig, self.min_val, self.max_val),
            interval=100,
            repeat_delay=1000,
            frames=len(range(0, dim[self.dimint])[:: self.step]),
        )
        if self.savegif is not None:
            ani.save(self.savegif, writer="pillow")
            self._show()
        else:
            self._show()

    def slicer(self):

        if self.dim in [None, "X", "x"]:
            plt.xlabel("Y")
            plt.ylabel("Z")
            plt.imshow(
                np.squeeze(self.cube[:, :, self.slice]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dim in ["Y", "y"]:
            plt.xlabel("X")
            plt.ylabel("Z")
            plt.imshow(
                np.squeeze(self.cube[:, self.slice]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dim in ["Z", "z"]:
            plt.xlabel("X")
            plt.ylabel("Y")
            plt.imshow(
                np.squeeze(self.cube[self.slice]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        self._show()

    def _show(self):
        if self.show_plot:
            plt.show()


plotimg = plotImg
