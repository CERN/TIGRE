import numpy as np
import re
import mrcfile

def readTltFile(fnTlt):
    try:
        with open(fnTlt, 'r') as tlt:
            lines = tlt.readlines()
            number_pattern = re.compile(r'^\s*([\d.-]+)')
            first_column = []
            for line in lines:
                match = number_pattern.match(line)
                if match:
                    first_column.append(float(match.group(1)))
            return np.array(first_column)
    except IOError:
        raise FileNotFoundError("Error parsing the tlt file")


def readTiltSeries(fnTs, normalizeTis):
    ts_aux = None#mrcfile.read(fnTs)
    with mrcfile.open(fnTs, permissive=True) as mrc:
        ts_aux = mrc.data.copy()

    dims = np.shape(ts_aux)
    xdim    = dims[2]
    ydim    = dims[1]
    nimages = dims[0]

    if normalizeTis == 'standard':
        for i in range(nimages):
            ti           = ts_aux[i, :, :]
            stdTi        = np.std(ti)
            meanTi       = np.mean(ti)
            ti           = (ti - meanTi) / stdTi
            ts_aux[i, :, :] = ti

    return ts_aux.astype(np.float32), xdim, ydim

def readXf(filename):
    """
    Parses an IMOD .xf file into a list of (matrix2x2, translation) pairs, one per line.

    Each line contains six floats: A11 A12 A21 A22 DX DY, where the 2x2 matrix
    controls rotation/scale/shear and (DX, DY) is the translation, both applied
    about the center of the image.
    """
    transforms = []
    with open(filename, 'r') as f:
        for lineNumber, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            values = line.split()
            if len(values) != 6:
                raise ValueError(
                    f"{filename}:{lineNumber}: expected 6 values, got {len(values)}"
                )

            a11, a12, a21, a22, dx, dy = (float(v) for v in values)
            matrix2x2 = np.array([[a11, a12], [a21, a22]])
            translation = np.array([dx, dy])
            transforms.append((matrix2x2, translation))

    return transforms
