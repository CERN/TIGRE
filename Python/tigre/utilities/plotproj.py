from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
def plotproj(projections):
    plt.ion()

    min_val = np.amin(projections)
    max_val = np.amax(projections)
    total_projections = projections.shape[0]
    for i in range(total_projections):
        plt.clf()
        plt.imshow(np.squeeze(projections[i]), cmap=plt.cm.gray, origin='lower', vmin=min_val, vmax=max_val)

        plt.gca().get_xaxis().set_ticks([])
        plt.gca().get_yaxis().set_ticks([])

        plt.xlabel("-> U")
        plt.ylabel("-> V")
        plt.title("Projection angle : " + str(i*360/total_projections))

        plt.colorbar()

        plt.pause(0.001)

def ppslice(projections,slice=None,Dim=2):
    if slice is None:
        slice = projections.shape[2]/2
    min_val = np.amin(projections)
    max_val = np.amax(projections)
    plt.clf()
    if Dim==0:
        plt.imshow(np.squeeze(projections[:,:,slice]), cmap=plt.cm.gray, origin='lower', vmin=min_val, vmax=max_val)
    if Dim==1:
        plt.imshow(np.squeeze(projections[:,slice]), cmap=plt.cm.gray, origin='lower', vmin=min_val, vmax=max_val)
    if Dim==2:
        plt.imshow(np.squeeze(projections[slice]), cmap=plt.cm.gray, origin='lower', vmin=min_val, vmax=max_val)
    plt.colorbar()
    plt.show()

