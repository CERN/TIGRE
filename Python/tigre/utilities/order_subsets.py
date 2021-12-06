from __future__ import division
import numpy as np


def _get_angularDistance_index(angles):
    # Sort values according to size
    angles_alpha = angles[:, 0]  # omit theta and psi
    index_alpha = np.argsort(angles_alpha)
    # Save the first values for index_alpha and angles in the new lists
    # (this needs to be done or argmin will return error otherwise)
    new_index = [index_alpha[0]]
    new_angles_alpha = [angles_alpha[index_alpha[0]]]
    # remove the first values from the original arrays
    angles_alpha = np.delete(angles_alpha[index_alpha], 0)
    index_alpha = np.delete(index_alpha, 0)
    # Run through list
    for i in range(len(angles_alpha)):
        # Run the comparison
        compangle = new_angles_alpha[i]
        angles_min = np.argmin(abs(abs(angles_alpha - compangle) - (np.pi / 2)))
        # Save the results to the new lists
        new_angles_alpha.append(angles_alpha[angles_min])
        new_index.append(index_alpha[angles_min])
        # delete old values from the original arrays
        angles_alpha = np.delete(angles_alpha, angles_min)
        index_alpha = np.delete(index_alpha, angles_min)

    return np.array(new_index)


def order_subsets(angles, blocksize, mode):
    """
    Returns
        block_alpha, index_alpha: lists of numpy.ndarray
    """

    # Parse no blocks
    if blocksize is None or blocksize == 1:
        index_alpha = np.arange(angles.shape[0])

        if mode is None or mode == "ordered":
            return list(angles), list(index_alpha)

        if mode == "random" or mode == "random2":
            np.random.shuffle(index_alpha)
            return list(angles[index_alpha]), list(index_alpha)
        if mode == "angularDistance":
            new_index = _get_angularDistance_index(angles)
            new_angles = angles[new_index]
            return list(np.array(new_angles, dtype=np.float32)), list(new_index)
        else:
            raise NameError("mode string not recognised: " + mode)

    # Parse with blocks
    elif blocksize > 1:
        nangles = angles.shape[0]
        block_count = (nangles+blocksize-1)//blocksize
        if mode == "random2":
            new_order = np.arange(nangles)
            np.random.shuffle(new_order)
            index_alpha = [new_order[i*blocksize:(i+1)*blocksize] for i in range(block_count-1)]
            index_alpha.append(new_order[(block_count-1)*blocksize:nangles])
            new_alpha = angles[new_order]
            block_alpha = [new_alpha[i*blocksize:(i+1)*blocksize] for i in range(block_count-1)]
            block_alpha.append(new_alpha[(block_count-1)*blocksize:nangles])
            return block_alpha, index_alpha

        if mode == "angularDistance":
            new_order = _get_angularDistance_index(angles)

            index_alpha = [new_order[i*blocksize:(i+1)*blocksize] for i in range(block_count-1)]
            index_alpha.append(new_order[(block_count-1)*blocksize:nangles])
            new_alpha = angles[new_order]
            block_alpha = [new_alpha[i*blocksize:(i+1)*blocksize] for i in range(block_count-1)]
            block_alpha.append(new_alpha[(block_count-1)*blocksize:nangles])
            return block_alpha, index_alpha

        # using list comprehension to form the blocks.
        oldindex = np.arange(nangles,dtype=np.int32)
        index_alpha = [oldindex[i*blocksize:(i+1)*blocksize] for i in range(0,block_count-1)]
        index_alpha.append(oldindex[(block_count-1)*blocksize:nangles])
        block_alpha = [angles[i*blocksize:(i+1)*blocksize] for i in range(0,block_count-1)]
        block_alpha.append(angles[(block_count-1)*blocksize:nangles])

        if mode is None or mode == "ordered":
            return block_alpha, index_alpha

        if mode == "random" or mode is None:
            new_order = np.arange(len(index_alpha))
            np.random.shuffle(new_order)

            index_alpha = [index_alpha[i] for i in new_order]
            block_alpha = [block_alpha[i] for i in new_order]
            return block_alpha, index_alpha
        else:
            raise NameError("mode string not recognised: " + mode)
