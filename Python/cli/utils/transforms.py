import numpy as np
import torch
import torch.nn.functional as F
from scipy import ndimage


def buildForwardMatricesXy(matrices2x2, translations, imageShape):
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
    GPU/batched equivalent of createMatrixFromXf + applyTransformation.
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

