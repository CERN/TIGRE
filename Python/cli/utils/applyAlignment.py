
import argparse

import numpy as np
import mrcfile
import torch
import torch.nn.functional as F
from io_utils import readXf, readTiltSeries
from scipy import ndimage

def createScipyMatrix(imageShape, angleDegrees, tx, ty, order="rotateFirst"):
    """
    Creates the 3x3 affine matrix and adapts it for scipy.ndimage (rows/columns + inverse matrix).
    """
    rows, columns = imageShape
    cx, cy = columns / 2.0, rows / 2.0

    theta = np.deg2rad(angleDegrees)
    cosT, sinT = np.cos(theta), np.sin(theta)

    # 1. Standard matrix in X,Y (same as before)
    R = np.array([
        [cosT, -sinT, cx * (1 - cosT) + cy * sinT],
        [sinT,  cosT, cy * (1 - cosT) - cx * sinT],
        [ 0.0,   0.0, 1.0]
    ])

    T = np.array([
        [1.0, 0.0, tx],
        [0.0, 1.0, ty],
        [0.0, 0.0, 1.0]
    ])

    if order == "rotateFirst":
        mXy = T @ R
    else:
        mXy = R @ T

    # 2. Conversion to SciPy coordinates (Row, Column) instead of (X, Y)
    # Swap rows and columns using a permutation matrix P
    P = np.array([
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0]
    ])
    mRowCol = P @ mXy @ P

    # 3. SciPy requires the INVERSE matrix for interpolation calculations
    mScipyInverse = np.linalg.inv(mRowCol)

    return mScipyInverse

def applyScipyTransformation(image, inverseMatrix3x3, fillValue=None):
    """
    Applies the transformation using scipy.ndimage.affine_transform.
    """
    if fillValue is None:
        fillValue = np.mean(image)

    # Extract the transformation sub-matrix (2x2) and the translation vector (offset)
    # SciPy expects the transformation matrix and the offset separately
    matrix2x2 = inverseMatrix3x3[:2, :2]
    offset = inverseMatrix3x3[:2, 2]

    # Apply the transformation
    # order=3 is cubic spline interpolation (good quality for cryo-EM)
    # order=5 would be quintic interpolation (higher quality, slower)
    transformedImage = ndimage.affine_transform(
        image,
        matrix=matrix2x2,
        offset=offset,
        order=3,
        mode='constant',
        cval=float(fillValue),
        prefilter=True # Important to keep True so splines don't blur the signal
    )

    return transformedImage

# def readMrc(filename):
#     """
#     Reads an MRC/MRCS file using mrcfile and returns its data as a numpy array.
#     """
#     with mrcfile.open(filename, permissive=True) as mrc:
#         data = mrc.data.copy()
#     return data


# def readXf(filename):
#     """
#     Parses an IMOD .xf file into a list of (matrix2x2, translation) pairs, one per line.

#     Each line contains six floats: A11 A12 A21 A22 DX DY, where the 2x2 matrix
#     controls rotation/scale/shear and (DX, DY) is the translation, both applied
#     about the center of the image.
#     """
#     transforms = []
#     with open(filename, 'r') as f:
#         for lineNumber, line in enumerate(f, start=1):
#             line = line.strip()
#             if not line:
#                 continue

#             values = line.split()
#             if len(values) != 6:
#                 raise ValueError(
#                     f"{filename}:{lineNumber}: expected 6 values, got {len(values)}"
#                 )

#             a11, a12, a21, a22, dx, dy = (float(v) for v in values)
#             matrix2x2 = np.array([[a11, a12], [a21, a22]])
#             translation = np.array([dx, dy])
#             transforms.append((matrix2x2, translation))

#     return transforms


def createScipyMatrixFromXf(imageShape, matrix2x2, translation):
    """
    Builds the 3x3 inverse matrix (scipy row/column convention) for an IMOD-style
    2x2 matrix and translation, both defined about the center of the image.
    """
    rows, columns = imageShape
    cx, cy = columns / 2.0, rows / 2.0

    a11, a12 = matrix2x2[0]
    a21, a22 = matrix2x2[1]
    dx, dy = translation

    mXy = np.array([
        [a11, a12, cx * (1 - a11) - a12 * cy + dx],
        [a21, a22, cy * (1 - a22) - a21 * cx + dy],
        [0.0, 0.0, 1.0]
    ])

    P = np.array([
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0]
    ])
    mRowCol = P @ mXy @ P

    return np.linalg.inv(mRowCol)


def buildForwardMatricesXy(matrices2x2, translations, imageShape):
    """
    Batched version of the (x, y) forward affine matrix built inside
    createScipyMatrixFromXf, one 3x3 matrix per transform. Uses the same
    center convention (columns / 2, rows / 2) as createScipyMatrixFromXf so
    both paths agree pixel-for-pixel.
    """
    rows, columns = imageShape
    cx, cy = columns / 2.0, rows / 2.0

    a = np.asarray(matrices2x2, dtype=np.float64)   # (N, 2, 2): [[a11, a12], [a21, a22]]
    t = np.asarray(translations, dtype=np.float64)  # (N, 2): [dx, dy]
    n = a.shape[0]

    m = np.zeros((n, 3, 3))
    m[:, :2, :2] = a
    m[:, 0, 2] = cx * (1 - a[:, 0, 0]) - a[:, 0, 1] * cy + t[:, 0]
    m[:, 1, 2] = cy * (1 - a[:, 1, 1]) - a[:, 1, 0] * cx + t[:, 1]
    m[:, 2, 2] = 1.0
    return m


def applyTorchTransforms(images, matrices2x2, translations, device="cuda",
                          interpolation="bicubic", batchSize=4):
    """
    GPU/batched equivalent of createScipyMatrixFromXf + applyScipyTransformation.
    Applies one (matrix2x2, translation) per image, batchSize images at a time
    in a single torch.nn.functional.grid_sample call, to bound GPU memory use
    on large stacks/images.

    Note: torch's "bicubic" mode is a cubic-convolution kernel, not the same
    cubic B-spline scipy's order=3 uses. Results are visually equivalent but
    not bit-identical to the scipy path.

    images: (N, H, W) array
    matrices2x2: (N, 2, 2) array
    translations: (N, 2) array
    """
    images = np.asarray(images)
    n, rows, columns = images.shape

    mForward = buildForwardMatricesXy(matrices2x2, translations, (rows, columns))
    mInverse = np.linalg.inv(mForward)  # output (x, y) -> input (x, y), pixel units

    # Output pixel grid, raw (x, y) coordinates, shared by every image
    yOut, xOut = torch.meshgrid(
        torch.arange(rows, dtype=torch.float32, device=device),
        torch.arange(columns, dtype=torch.float32, device=device),
        indexing="ij"
    )
    pOut = torch.stack([xOut, yOut, torch.ones_like(xOut)], dim=0).reshape(3, -1)  # (3, H*W)

    output = np.empty_like(images)
    with torch.no_grad():
        for start in range(0, n, batchSize):
            end = min(start + batchSize, n)
            batchImages = torch.as_tensor(images[start:end], dtype=torch.float32, device=device)
            batchInverseT = torch.as_tensor(mInverse[start:end], dtype=torch.float32, device=device)

            # Map every output pixel to its input pixel, per image in the batch
            pIn = torch.einsum("nij,jk->nik", batchInverseT, pOut)  # (B, 3, H*W)
            xIn, yIn = pIn[:, 0], pIn[:, 1]

            # Normalization to [-1, 1] matching grid_sample's align_corners=False
            # convention (pixel i spans [i, i+1), norm = (2*pixel + 1) / size - 1),
            # so it lines up with the columns/2, rows/2 pixel-center convention above.
            b = end - start
            xInNorm = ((2.0 * xIn + 1.0) / columns - 1.0).reshape(b, rows, columns)
            yInNorm = ((2.0 * yIn + 1.0) / rows - 1.0).reshape(b, rows, columns)
            grid = torch.stack([xInNorm, yInNorm], dim=-1)

            batchOutput = F.grid_sample(
                batchImages.unsqueeze(1), grid, mode=interpolation,
                padding_mode="zeros", align_corners=False
            ).squeeze(1)

            # Emulate scipy's mode='constant', cval=mean(image) for out-of-bounds samples
            fillValues = batchImages.mean(dim=(1, 2)).view(b, 1, 1)
            outOfBounds = (grid[..., 0].abs() > 1) | (grid[..., 1].abs() > 1)
            batchOutput = torch.where(outOfBounds, fillValues.expand_as(batchOutput), batchOutput)

            output[start:end] = batchOutput.cpu().numpy()

    return output


def parseArgs():
    parser = argparse.ArgumentParser(
        description="Apply a stack of IMOD .xf 2D transformations to an MRC image stack."
    )
    parser.add_argument("-i", dest="input", required=True,
                         help="Input MRC/MRCS image stack")
    parser.add_argument("-xf", dest="xf", required=True,
                         help="IMOD .xf file with one transformation per image")
    parser.add_argument("-o", dest="output", required=True,
                         help="Output interpolated MRCS image stack")
    parser.add_argument("--device", dest="device", default=None,
                         choices=["cuda", "cpu"],
                         help="Device for the interpolation (default: cuda if available, else cpu)")
    return parser.parse_args()


def main():
    args = parseArgs()
    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")

    images = readTiltSeries(args.input, normalizeTis = False)
    if images.ndim == 2:
        images = images[np.newaxis, ...]

    transforms = readXf(args.xf)
    if len(transforms) != images.shape[0]:
        raise ValueError(
            f"Number of transformations ({len(transforms)}) does not match "
            f"number of images ({images.shape[0]})"
        )

    matrices2x2 = np.stack([m for m, _ in transforms])
    translations = np.stack([t for _, t in transforms])

    output = applyTorchTransforms(images, matrices2x2, translations, device=device)
    output = output.astype(np.float32)

    with mrcfile.new(args.output, overwrite=True) as mrc:
        mrc.set_data(output)


if __name__ == "__main__":
    main()
