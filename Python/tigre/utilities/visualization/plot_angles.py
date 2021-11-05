""" Visualization of ordered subset of angles
#   This file is part of the TIGRE Toolbox

#   Copyright (c) 2015, University of Bath and
#                       CERN-European Organization for Nuclear Research
#                       All rights reserved.

#   License:            Open Source under BSD.
#                       See the full license at
#                       https://github.com/CERN/TIGRE/license.txt

#   Contact:            tigre.toolbox@gmail.com
#   Codes:              https://github.com/CERN/TIGRE/
# --------------------------------------------------------------------------
#   Coded by:          Tomoyuki SADAKANE
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

class PlotAngles:
    """Generate ordered subset angles animation
    Parameters
    ------------
        angle_blocks: List of subsets of extended (alpha, theta, psi) angles.
        angle_index: List of subsets of indeces of angles.
        savegif: Specify file name to create a gif file.
        show_plot: Set False not to display the plot.
        figsize: figsize to be passed to plt.figure()
        interval: Animation interval of the frames in ms.
        draw_number: Draw number in the plot. Default is True.
        size_marker: Size of the angle marker. Default is None to use the default of plt.scatter().
        angles_extended: If angles are not extended, set False. Default is True.
    """
    def __init__(self,
                angle_blocks,
                angle_index,
                savegif=None,
                show_plot=True,
                figsize=None,
                interval=500,
                draw_number=True,
                size_marker=None,
                angles_extended=True):
        if not angles_extended:
            if type(angle_blocks) is list:
                for idx, block in enumerate(angle_blocks):
                    a1 = block
                    a2 = np.zeros(a1.shape[0], dtype=np.float32)
                    angle_blocks[idx] = np.vstack((a1, a2, a2)).T
            elif type(angle_blocks) is np.ndarray:
                if angle_blocks.ndim == 1:
                    # treat as blocksize = 1
                    a1 = angle_blocks
                    a2 = np.zeros(a1.shape[0], dtype=np.float32)
                    angle_blocks = np.vstack((a1, a2, a2)).T
                    angle_blocks = list(angle_blocks)

        if type(angle_blocks) is np.ndarray:
            angle_blocks = [angle_blocks]
        if type(angle_index) is np.ndarray:
            angle_index = [angle_index]
        if len(angle_blocks) == 0 or len(angle_index) == 0:
            raise TypeError("angle_blocks and angle_index must be a lists of sizes greater than 0.")

        if angle_blocks[0].ndim != angle_index[0].ndim:
            # angles is extended
            if angle_blocks[0].ndim == 1:
                self.angle_blocks = [np.array([angle_blocks[i][0]]) for i in range(len(angle_blocks))]
            else:
                self.angle_blocks = [angle_blocks[i][:, 0] for i in range(len(angle_blocks))]
        else:
            self.angle_blocks = angle_blocks

        self.angle_index = angle_index

        self.show_plot  = show_plot
        self.figsize = figsize
        self.savegif = savegif

        self.interval = interval
        self.draw_number = draw_number
        self.size_marker = size_marker

        self._run()

    def _update_frame(self, it, fig, radius):
        fig.clf()
        axis = fig.add_subplot(1, 1, 1)
        rangewidth = radius*1.1
        axis.set_aspect("equal")
        axis.set_xlim([-rangewidth, +rangewidth]*radius)
        axis.set_ylim([-rangewidth, +rangewidth]*radius)
        axis.axhline(0, ls="solid", color="gray", linewidth = 0.5)
        axis.axvline(0, ls="solid", color="gray", linewidth = 0.5)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        radius_text = radius *0.9
        block_count = len(self.angle_blocks)
        for idx in range(it):
            if idx < block_count:
                for idx_block in range(self.angle_blocks[idx].shape[0]):
                    x = radius*np.cos(self.angle_blocks[idx][idx_block])
                    y = radius*np.sin(self.angle_blocks[idx][idx_block])
                    axis.scatter(x, y, color="red", s=self.size_marker)
                    if self.draw_number:
                        x = radius_text*np.cos(self.angle_blocks[idx][idx_block])
                        y = radius_text*np.sin(self.angle_blocks[idx][idx_block])
                        axis.text(x, y, str(idx), verticalalignment = "center", horizontalalignment = "center")

        if it < block_count:
            for idx_block in range(self.angle_blocks[it].shape[0]):
                self.x_old = radius*np.cos(self.angle_blocks[it][idx_block])
                self.y_old = radius*np.sin(self.angle_blocks[it][idx_block])
                axis.scatter(self.x_old, self.y_old, color="green", s=self.size_marker)
                if self.draw_number:
                    x = radius_text*np.cos(self.angle_blocks[it][idx_block])
                    y = radius_text*np.sin(self.angle_blocks[it][idx_block])
                    axis.text(x, y, str(it), verticalalignment = "center", horizontalalignment = "center")


    def _run(self):

        radius = 1
        self.x_old = radius*np.cos(self.angle_blocks[0][0])
        self.y_old = radius*np.sin(self.angle_blocks[0][0])

        fig = plt.figure(figsize=self.figsize)
        ani = animation.FuncAnimation(
                            fig,
                            self._update_frame,
                            fargs=(fig, radius),
                            interval=self.interval,
                            repeat_delay=1000,
                            frames = len(self.angle_blocks)+1)
        if self.savegif is not None:
            ani.save(self.savegif, writer="pillow")
            self._show()
        else:
            self._show()

    def _show(self):
        if self.show_plot:
            plt.show()

plot_angles = PlotAngles
