import numpy as np
import pytest
from tigre.utilities.im_3d_denoise import im3ddenoise


@pytest.mark.parametrize("shape", [(256, 256), (256,)])
def test_non_3d_input_raises_value_error(shape):
    """im3ddenoise needs a 3D array and says so (issue #681)."""
    img = np.random.rand(*shape).astype(np.float32)
    with pytest.raises(ValueError, match="3D array"):
        im3ddenoise(img, iter=5)
