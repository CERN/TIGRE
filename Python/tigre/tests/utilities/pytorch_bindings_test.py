from typing import Union
import tigre
import torch
import pytest 
import numpy as np

import tigre.utilities.gpu as gpu
from tigre.utilities.geometry import Geometry
from temp_pytorch_bindings import A, At

# We begin by defining the pytest fixtures
# see https://docs.pytest.org/en/6.2.x/fixture.html
    
projectors = []
projectors.extend(
    (pytest.param(proj_cfg)
     for proj_cfg in [
        'cone2D uniform',
        'cone3D uniform'
        ])
)

projector_ids = [
    " geometry_name='{}' angles='{}' ".format(*p.values[0].split())
    for p in projectors
]

@pytest.fixture(scope='function', params=projectors, ids=projector_ids)
def projector(request) -> Union[Geometry, np.ndarray]:
    geometry_name, angles = request.param.split()
    if angles == 'uniform':
        angles = np.linspace(0, np.pi, 200)
    else:
        raise ValueError('angles not valid')
    
    if geometry_name == 'cone2D':
        volume   = np.ones((1, 256, 256), dtype=np.float32)
        geo = tigre.geometry()
        geo.DSD = 1536  # Distance Source Detector      (mm)
        geo.DSO = 1000  # Distance Source Origin        (mm)
        # Image parameters
        geo.nVoxel = np.array([1, 256, 256])  # number of voxels              (vx)
        geo.sVoxel = np.array([1, 256, 256])  # total size of the image       (mm)
        geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
        # Detector parameters
        geo.nDetector = np.array([1, 512])  # number of pixels              (px)
        geo.dDetector = np.array([geo.dVoxel[0], 0.8])  # size of each pixel            (mm)
        geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm)
        # Offsets
        geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
        geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
        # MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!
        geo.mode = "cone"
        
    elif geometry_name == 'cone3D':
        volume   = np.ones((64, 64, 64), dtype=np.float32)
        geo = tigre.geometry()
        geo.DSD = 1536  # Distance Source Detector      (mm)
        geo.DSO = 1000  # Distance Source Origin        (mm)
        # Image parameters
        geo.nVoxel = np.array([64, 64, 64])  # number of voxels              (vx)
        geo.sVoxel = np.array([64, 64, 64])  # total size of the image       (mm)
        geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
        # Detector parameters
        geo.nDetector = np.array([32, 32])  # number of pixels              (px)
        geo.dDetector = np.array([geo.dVoxel[0], 0.8])  # size of each pixel            (mm)
        geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm)
        # Offsets
        geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
        geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
        # MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!
        geo.mode = "cone"
    else:
        raise ValueError('geometry_name not valid')
    
    geo.angles = angles

    return geo, volume

@pytest.fixture(params=gpu.getGpuNames())
def tigre_device(projector, request):
    geometry = projector[0]
    n_vox = geometry.nVoxel[0]
    if n_vox == 1:
        return gpu.getGpuIds(request.param)[0]
    else:
        return gpu.getGpuIds(request.param)
 
@pytest.fixture(params=[i for i in range(len(gpu.getGpuNames()))])
def torch_device(request):
    return request.param

def test_angles_assertion(projector, tigre_device):
    geometry, numpy_volume = projector
    geometry.angles = None
    with pytest.raises(AssertionError, match='Initialise the angles'):
        forward_operator = A(geometry, tigre_device)
        raise ValueError('Initialise the angles')

def test_forward_operator(projector, tigre_device, torch_device):
    geometry, numpy_volume = projector

    torch_volume = torch.from_numpy(numpy_volume).unsqueeze(0).unsqueeze(0).to(f'cuda:{torch_device}')
    numpy_sinogram = tigre.Ax(
        numpy_volume, 
        geometry, 
        geometry.angles, 
        gpuids=tigre_device
        )
    
    forward_operator = A(geometry, tigre_device)
    torch_sinogram = forward_operator(torch_volume)[0,0].detach().cpu().numpy()
    np.testing.assert_array_almost_equal_nulp(numpy_sinogram, torch_sinogram)

def test_backward_operator(projector, tigre_device, torch_device):
    geometry, numpy_volume = projector    
    numpy_sinogram = tigre.Ax(
        numpy_volume, 
        geometry, 
        geometry.angles, 
        gpuids=tigre_device
        )
    numpy_volume = tigre.Atb(
        numpy_sinogram, 
        geometry, 
        geometry.angles, 
        gpuids=tigre_device
        )
    torch_sinogram = torch.from_numpy(numpy_sinogram).unsqueeze(0).unsqueeze(0).to(f'cuda:{torch_device}')
    
    backward_operator = At(geometry, tigre_device)
    torch_volume = backward_operator(torch_sinogram)[0,0].detach().cpu().numpy()
    np.testing.assert_array_almost_equal_nulp(numpy_volume, torch_volume)

def test_forward_operator_requires_grad(projector, tigre_device, torch_device):
    geometry, numpy_volume = projector

    torch_volume = torch.from_numpy(numpy_volume).unsqueeze(0).unsqueeze(0).to(f'cuda:{torch_device}').requires_grad_(True)
    
    forward_operator = A(geometry, tigre_device)

    sinogram:torch.Tensor = forward_operator(torch_volume)

    assert sinogram.requires_grad

def test_backward_operator_requires_grad(projector, tigre_device, torch_device):
    geometry, numpy_volume = projector

    numpy_sinogram = tigre.Ax(
        numpy_volume, 
        geometry, 
        geometry.angles, 
        gpuids=tigre_device
        )
    torch_sinogram = torch.from_numpy(numpy_sinogram).unsqueeze(0).unsqueeze(0).to(f'cuda:{torch_device}').requires_grad_(True)
    
    backward_operator = At(geometry, tigre_device)

    volume:torch.Tensor = backward_operator(torch_sinogram)

    assert volume.requires_grad

def test_gradient_forward_operator(projector, tigre_device, torch_device):
    geometry, numpy_volume = projector
    torch_volume = torch.from_numpy(numpy_volume).unsqueeze(0).to(f'cuda:{torch_device}')
    sinogram = tigre.Ax(
        numpy_volume, 
        geometry, 
        geometry.angles,
        gpuids = tigre_device
        )
    sinogram = torch.from_numpy(sinogram).to(f'cuda:{torch_device}')
    torch_volume.requires_grad_(True)

    forward_operator = A(geometry, tigre_device)

    torch_sinogram = forward_operator(torch_volume)

    loss_f = torch.nn.MSELoss()

    L = loss_f(torch_sinogram[0], sinogram)

    L.backward()

    assert torch_volume.grad.size() == torch_volume.size()


def test_gradient_backward_operator(projector, tigre_device, torch_device):
    geometry, numpy_volume = projector
    sinogram = tigre.Ax(
        numpy_volume, 
        geometry, 
        geometry.angles,
        gpuids = tigre_device
        )
    volume = tigre.Atb(
        sinogram, 
        geometry, 
        geometry.angles,
        gpuids = tigre_device
        )
    
    volume = torch.from_numpy(volume).to(f'cuda:{torch_device}')
    torch_sinogram = torch.from_numpy(sinogram).to(f'cuda:{torch_device}').unsqueeze(0)
    torch_sinogram.requires_grad_(True)

    backward_operator = At(geometry, tigre_device)

    torch_volume = backward_operator(torch_sinogram)

    loss_f = torch.nn.MSELoss()

    L = loss_f(torch_volume[0], volume)

    L.backward()

    assert torch_sinogram.grad.size() == torch_sinogram.size()