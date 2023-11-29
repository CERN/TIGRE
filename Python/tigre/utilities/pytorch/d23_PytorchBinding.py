import numpy as np

import tigre
import Operator
from tigre.utilities import sl3d, gpu

import torch

import matplotlib.pyplot as plt


def get_default_geometry():
    geo = tigre.geometry()
    #%% Geometry
    geo.DSD = 1536  # Distance Source Detector      (mm) #type:ignore
    geo.DSO = 1000  # Distance Source Origin        (mm) #type:ignore
    # Image parameters
    geo.nVoxel = np.array([1, 256, 256])  # number of voxels              (vx)
    geo.sVoxel = np.array([1, 256, 256])  # total size of the image       (mm)
    geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
    # Detector parameters
    geo.nDetector = np.array([1, 512])  # number of pixels              (px)
    geo.dDetector = np.array([geo.dVoxel[0], 0.8])  # size of each pixel            (mm)
    geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm) #type:ignore
    # Offsets
    geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
    geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
    # MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!
    geo.mode = "parallel"
    return geo

if __name__ == '__main__':
    plt.style.use('grayscale')


    #%% Define angles of projection and load phantom image
    # define projection angles (in radians)
    angles = np.linspace(0, 2 * np.pi, 50)
     #%% Set GPU ids (we assume only one GPU is available)
    gpu_index = 2
    gpuids = gpu.getGpuIds()
    gpuids.devices = [gpu_index]

    two_d = False
    if two_d:
        geo = get_default_geometry()
        x = torch.zeros((1,1,256,256), device='cuda:2', requires_grad=True)
        x_target = torch.from_numpy(np.load('test_phantom_2D.npy').astype(np.float32)[::2,::2])
        x_target = x_target.to('cuda:2').unsqueeze(0).unsqueeze(0)
        plt.matshow(x_target[0,0].detach().cpu())
        plt.savefig(f'x_target.jpg')
        plt.clf()
    else:
        #%% Define default geometry
        geo = tigre.geometry_default(high_resolution=False)
        # load phantom image
        phantom_type = (
            "yu-ye-wang"  # Default of Python TIGRE Shepp-Logan phantom. Improved visual perception
        )
        shepp_phantom:np.ndarray = sl3d.shepp_logan_3d(
            geo.nVoxel, phantom_type=phantom_type
        )
        print(f'Phantom shape: {shepp_phantom.shape}')
        x = torch.zeros(
        size = (1,1)+shepp_phantom.shape,
        device=f'cuda:{gpu_index}',
        requires_grad=True
        )
        x_target = torch.from_numpy(shepp_phantom).to(f'cuda:{gpu_index}').unsqueeze(0).unsqueeze(0)
        plt.matshow(x_target[0,0,32].detach().cpu())
        plt.savefig(f'x_target.jpg')
        plt.clf()

    #%% Define torch bindings
    # Perform a L2-loss minimisation of the data discrepancy term
    loss = torch.nn.MSELoss()
    assert x.size() == x_target.size()
    # Define the operator
    operator = Operator.OperatorModule(geo, angles, gpuids)

    # Define the optimiser
    optimiser = torch.optim.SGD([x], lr=1)

    n_steps = 201

    '''sinogram = operator(x_target)
    if two_d:
        plt.matshow(sinogram[0,0].detach().cpu())
        plt.savefig(f'measurement.jpg')
        plt.clf()
    else:
        plt.matshow(sinogram[0,0,0].detach().cpu())
        plt.savefig(f'measurement.jpg')
        plt.clf()'''

    for i in range(n_steps):
        optimiser.zero_grad()
        mse_loss = loss(operator(x), operator(x_target))
        mse_loss.backward()
        optimiser.step()
        print(i, mse_loss.item())

    if two_d:
        plt.matshow(x[0,0].detach().cpu())
        plt.savefig(f'reconstruction.jpg')
        plt.clf()
    else:
        plt.matshow(x[0,0,0].detach().cpu())
        plt.savefig(f'reconstruction.jpg')
        plt.clf()



