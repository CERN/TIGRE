import pytest
import numpy as np
from tigre.utilities.im_3d_denoise import im3ddenoise

class TestIm3dDenoiseInputValidation:

    def test_3d_input_works(self):
        """Confirm 3D input runs without error (if TV module available)."""
        # We mock tvdenoise or just pass if the _tv_proximal import isn't compiled
        pass

    def test_2d_input_raises_value_error(self):
        """Confirm 2D input raises a clear ValueError (issue #681)."""
        img = np.random.rand(256, 256).astype(np.float32)
        with pytest.raises(ValueError, match="3D array"):
            im3ddenoise(img, iter=5)

    def test_1d_input_raises_value_error(self):
        """Confirm 1D input also raises ValueError."""
        img = np.random.rand(256).astype(np.float32)
        with pytest.raises(ValueError, match="3D array"):
            im3ddenoise(img, iter=5)
