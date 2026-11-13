import numpy as np
import pytest
from tigre.design_api import linear_source_surface

def test_linear_source_exact_endpoints_and_fixed_detector():
    pixels=np.zeros((2,3,3)); pixels[...,0]=-10
    surface=linear_source_surface(pixels,[10,-2,0],[10,2,0],5)
    np.testing.assert_allclose(surface.source_positions_mm,np.column_stack((np.full(5,10),np.linspace(-2,2,5),np.zeros(5))))
    np.testing.assert_array_equal(surface.pixel_centers_mm,np.broadcast_to(pixels,(5,2,3,3)))
    with pytest.raises(ValueError):
        linear_source_surface(pixels,[10,-2,0],[10,2,0],0)
