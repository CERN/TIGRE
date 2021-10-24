from __future__ import print_function

import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


class plotProj:
    # NOTE: type help(plotImg) after importing in order to get a readable manual.
    (
        "\n"
        "plotProj(proj, dim) \n"
        "    plots figure \n"
        "default: progressive in slices following\n"
        "    axis (dim)\n"
        "Parameters \n"
        "---------- \n"
        "proj : Any 3D numpy array \n"
        "\n"
        'dim : ("U","V","T","u","v","t"), optional \n'
        '       default is "T"\n'
        "       NOTE: string arguments!"
        "\n"
        "angles: Any 1D numpy array. \n"
        "        Its length must be the same as proj.shape[0].\n"
        '        Works only when dim is "T" or "t"\n'
        "slice: int, optional\n"
        "     returns page of matrix according to index\n"
        "step: int, optional\n"
        "      Sets the step size between slice and slice."
        "      Step is 1 by default.\n"
        "savegif: string, optional\n"
        "         Saves the image as .gif with the file name\n"
        "show_plot: bool, optional\n"
        "           Sets whether to show the plot.\n"
        "           Default is None, automatically detects matplotlib backend\n"
        "           and decides whether to call plt.show.\n"
        "Examples:\n"
        "---------\n"
        "a=np.ones([3,3,3])\n"
        "plotImg(a)\n"
        ">>>returns plot along dim T\n"
        'plotImg(a,dim="v")\n'
        ">>>returns plot along dim V\n"
    )

    def __init__(
        self,
        proj,
        angles=None,
        dim=None,
        slice=None,
        step=1,
        savegif=None,
        colormap="gray",
        clims=None,
        show_plot=None,
    ):
        self.proj = proj
        self.dim = dim
        self.slice = slice
        self.dimint = None  # keeps track of what dim
        self.dimlist = ["U", "V", "T", "u", "v", "t", None]  # accepted parameters for dim
        self.step = step
        self.savegif = savegif
        self.angles = angles
        self.colormap = colormap
        if clims is None:
            self.min_val = np.amin(self.proj)
            self.max_val = np.amax(self.proj)
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
        if self.step is None or self.step == 0:
            self.step = 1
        if self.savegif == "":
            self.savegif == None
        if self.slice is None:
            self.run()
        if self.slice is not None:
            self.slicer()

    def run(self):
        if self.dim not in self.dimlist and self.dim is not None:
            raise NameError("check inputs for dim, should be string.")
        if self.angles is not None and self.angles.shape[0] != self.proj.shape[0]:
            raise NameError("check inputs for angles, should be size of proj.shape[0]")

        if self.dim in ["U", "u"]:
            self.dimint = 2
            self.dimlist = ["->V", "->T", "U"]
            self.run_plot()
        if self.dim in ["V", "v"]:
            self.dimint = 1
            self.dimlist = ["->U", "->T", "V"]
            self.run_plot()
        if self.dim in [None, "T", "a"]:
            self.dimint = 0
            self.dimlist = ["->U", "->V", "T"]
            self.run_plot()

    def update_frame(self, it, fig, min_val, max_val):
        i = range(0, self.proj.shape[self.dimint])[:: self.step][it]
        fig.clf()
        axis = fig.add_subplot(1, 1, 1)
        if self.dimint == 2:
            mappable = axis.imshow(
                np.squeeze(self.proj[:, :, i]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dimint == 1:
            mappable = axis.imshow(
                np.squeeze(self.proj[:, i]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dimint == 0:
            mappable = axis.imshow(
                np.squeeze(self.proj[i]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        # axis.get_xaxis().set_ticks([])
        # axis.get_yaxis().set_ticks([])
        axis.set_xlabel(self.dimlist[0])
        axis.set_ylabel(self.dimlist[1])
        if self.angles is not None:
            axis.set_title(
                "{}:{}, alpha={:+.3f} pi".format(self.dimlist[2], i, self.angles[i] / np.pi)
            )
        else:
            axis.set_title(self.dimlist[2] + ":" + str(i))
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(mappable, cax=cax)
        # plt.pause(0.01)

    def run_plot(self):

        dim = self.proj.shape

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

        if self.dim in ["U", "u"]:
            plt.xlabel("V")
            plt.ylabel("T")
            plt.imshow(
                np.squeeze(self.proj[:, :, self.slice]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dim in ["V", "v"]:
            plt.xlabel("U")
            plt.ylabel("T")
            plt.imshow(
                np.squeeze(self.proj[:, self.slice]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        if self.dim in [None, "T", "t"]:
            if self.angles is not None:
                plt.title("alpha={:+.3f} pi".format(self.angles[self.slice] / np.pi))
            plt.xlabel("U")
            plt.ylabel("V")
            plt.imshow(
                np.squeeze(self.proj[self.slice]),
                cmap=self.colormap,
                origin="lower",
                vmin=self.min_val,
                vmax=self.max_val,
            )
        self._show()

    def _show(self):
        if self.show_plot:
            plt.show()


def plotSinogram(proj, posV, show_plot=None):  # noqa: N803
    """
    plotSinogram(proj, posV)
        plots sinogram at V=posV

    Parameters
    ----------
    proj : Any 3D numpy array
    posV : integer. in range of 0:proj.shape[1].
    """
    plotProj(proj, dim="V", slice=posV, show_plot=show_plot)


plotproj = plotProj
