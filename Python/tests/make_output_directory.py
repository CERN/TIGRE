import os

import numpy as np


if __name__ == "__main__":
    dirname = os.path.dirname(__file__)
    if "logs" not in os.listdir(dirname):
        os.system("mkdir " + os.path.join(dirname, "logs"))
    outdir = os.path.join(dirname, "logs")

    outdir = os.path.join(outdir, "test" + str(len(os.listdir(outdir))))

    os.system("mkdir " + outdir)

    np.save(os.path.join(dirname, "targetdir.npy"), outdir)
