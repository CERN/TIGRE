from __future__ import print_function

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plotproj(projections):
    plt.ion()

    min_val = np.amin(projections)
    max_val = np.amax(projections)
    total_projections = projections.shape[0]
    for i in range(total_projections):
        plt.clf()
        plt.imshow(
            np.squeeze(projections[i]), cmap=plt.cm.gray, origin="lower", vmin=min_val, vmax=max_val
        )

        plt.gca().get_xaxis().set_ticks([])
        plt.gca().get_yaxis().set_ticks([])

        plt.xlabel("-> U")
        plt.ylabel("-> V")
        plt.title("Projection angle : " + str(i * 360 / total_projections))

        plt.colorbar()

        plt.pause(0.001)


def ppslice(projections, slice=None, Dim=2):  # noqa: N803
    if slice is None:
        slice = projections.shape[2] / 2
    min_val = np.amin(projections)
    max_val = np.amax(projections)
    plt.clf()
    if Dim == 0:
        plt.imshow(
            np.squeeze(projections[:, :, slice]),
            cmap=plt.cm.gray,
            origin="lower",
            vmin=min_val,
            vmax=max_val,
        )
    if Dim == 1:
        plt.imshow(
            np.squeeze(projections[:, slice]),
            cmap=plt.cm.gray,
            origin="lower",
            vmin=min_val,
            vmax=max_val,
        )
    if Dim == 2:
        plt.imshow(
            np.squeeze(projections[slice]),
            cmap=plt.cm.gray,
            origin="lower",
            vmin=min_val,
            vmax=max_val,
        )
    plt.colorbar()
    plt.show()


class plotProj:  # noqa: N801
    """
    plotProj(proj, dim)
        plots figure
    default: progressive in slices following
        axis (dim)

    Parameters
    ----------
    proj : Any 3D numpy array

    dim : ("U","V","T","u","v","t"), optional
        default is "T"
        NOTE: string arguments!

    angles: Any 1D numpy array.
            Its length must be the same as proj.shape[0].
            Works only when dim is "T" or "t"

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
    >>>returns plot along dim T

    plotImg(a,dim="v")
    >>>returns plot along dim V
    """

    def __init__(self, proj, angles=None, dim=None, slice=None, step=1, savegif=None):
        self.proj = proj
        self.dim = dim
        self.slice = slice
        self.dimint = None  # keeps track of what dim
        self.dimlist = ["U", "V", "T", "u", "v", "t", None]  # accepted parameters for dim
        self.step = step
        self.savegif = savegif
        self.angles = angles
        if self.step is None or self.step == 0:
            self.step = 1
        if self.savegif == "":
            self.savegif = None
        if self.slice is None:
            self.run()
        if self.slice is not None:
            self.slicer()

    def run(self):
        if self.dim not in self.dimlist:
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
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dimint == 1:
            mappable = axis.imshow(
                np.squeeze(self.proj[:, i]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dimint == 0:
            mappable = axis.imshow(
                np.squeeze(self.proj[i]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
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
        min_val = np.amin(self.proj)
        max_val = np.amax(self.proj)
        dim = self.proj.shape

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
        min_val = np.amin(self.proj)
        max_val = np.amax(self.proj)
        if self.dim in ["U", "u"]:
            plt.xlabel("V")
            plt.ylabel("T")
            plt.imshow(
                np.squeeze(self.proj[:, :, self.slice]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dim in ["V", "v"]:
            plt.xlabel("U")
            plt.ylabel("T")
            plt.imshow(
                np.squeeze(self.proj[:, self.slice]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        if self.dim in [None, "T", "t"]:
            if self.angles is not None:
                plt.title("alpha={:+.3f} pi".format(self.angles[self.slice] / np.pi))
            plt.xlabel("U")
            plt.ylabel("V")
            plt.imshow(
                np.squeeze(self.proj[self.slice]),
                cmap=plt.cm.gray,
                origin="lower",
                vmin=min_val,
                vmax=max_val,
            )
        plt.show()


def plotSinogram(proj, posV):  # noqa: N803
    """
    plotSinogram(proj, posV)
        plots sinogram at V=posV

    Parameters
    ----------
    proj : Any 3D numpy array
    posV : integer. in range of 0:proj.shape[1].
    """
    plotProj(proj, dim="V", slice=posV)
