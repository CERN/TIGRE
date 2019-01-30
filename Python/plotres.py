from matplotlib import pyplot as plt
import numpy as np
import os
import sys


for filename in os.listdir(os.curdir):
    if filename.endswith('.npy'):
        mat = np.load(filename)
        plt.matshow(mat[mat.shape[0]/2])
        plt.title(filename)
        plt.colorbar()
plt.show()
