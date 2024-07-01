# DEMO 23: ### Running TIGRE's PyTorch bindings 
#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# This file is part of the TIGRE Toolbox
#
# Copyright (c) 2015, University of Bath and
#                     CERN-European Organization for Nuclear Research
#                     All rights reserved.
#
# License:            Open Source under BSD.
#                     See the full license at
#                     https://github.com/CERN/TIGRE/blob/master/LICENSE
#
# Contact:            tigre.toolbox@gmail.com
# Codes:              https://github.com/CERN/TIGRE/
# Coded by:           Emilien Valat
# --------------------------------------------------------------------------
#
import tigre
import torch

import numpy as np
import tigre.utilities.gpu as gpu
from tigre.utilities.pytorch_bindings import A, At

def get_default_3Dgeometry():
    geo = tigre.geometry()
    # Geometry
    geo.DSD = 1536  # Distance Source Detector      (mm) #type:ignore
    geo.DSO = 1000  # Distance Source Origin        (mm) #type:ignore
    # Image parameters
    geo.nVoxel = np.array([256, 256, 256])  # number of voxels              (vx)
    geo.sVoxel = np.array([256, 256, 256])  # total size of the image       (mm)
    geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
    # Detector parameters
    geo.nDetector = np.array([256,256])  # number of pixels              (px)
    geo.dDetector = np.array([geo.dVoxel[0], 1])  # size of each pixel            (mm)
    geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm) #type:ignore
    # Offsets
    geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
    geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
    # MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!
    geo.mode = "parallel"
    return geo

def get_default_2Dgeometry():
    geo = tigre.geometry()
    # Geometry
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

def test_2D_operators(
        n_iterations = 5
    ):
    ### 2D test
    print('Running 2D test: minimising data discrepancy')
    ### Instanciate the volumes
    print('Input  Volume: Tensor of size [2,1,256,256]')
    print('Target Volume: Tensor of size [2,1,256,256]')
    input_volume  = torch.nn.Parameter(torch.zeros((2,1,256,256), device=f'cuda:{PYTORCH_GPU_ID}', requires_grad=True))
    target_volume = torch.ones((2,1,256,256), device=f'cuda:{PYTORCH_GPU_ID}', requires_grad=True)
    
    ### Instanciate the geometry
    geo = get_default_2Dgeometry()

    ### Instanciate the differentiable modules
    forward_operator = A(
        geo, np.linspace(0, np.pi, 200), TIGRE_GPU_ID
    )
    backward_operator = At(
        geo, np.linspace(0, np.pi, 200), TIGRE_GPU_ID
    )

    ### Instanciate a simple convolutional layer, the loss and the optimiser
    conv = torch.nn.Conv2d(1,1,1).to(device=torch.device(f'cuda:{PYTORCH_GPU_ID}'))
    mse = torch.nn.MSELoss()
    optimiser = torch.optim.Adam(
        [input_volume] + list(conv.parameters()), 1e-4
    )
    
    ### Compute sinograms
    input_sinogram  = forward_operator(input_volume)
    target_sinogram = forward_operator(target_volume)

    print(f'Testing the forward operator Ax for {n_iterations} iterations')
    for i in range(10):
        optimiser.zero_grad()
        approximation = conv(forward_operator(input_volume))
        loss = mse(approximation, target_sinogram)
        print(f'MSE Loss: {loss.item():.3f}')
        loss.backward()
        optimiser.step()
    print('Forward operator tests done!')

    print(f'Testing the backward operator Atb for {n_iterations} iterations')
    for i in range(10):
        optimiser.zero_grad()
        approximation = conv(backward_operator(input_sinogram))
        loss = mse(approximation, target_volume)
        print(f'MSE Loss: {loss.item():.5f}')
        loss.backward()
        optimiser.step()
    print('Backward operator tests done!')

def test_3D_operators(
        n_iterations = 5
    ):
    ### 3D test
    print('Running 3D test: minimising data discrepancy')
    ### Instanciate the volumes
    print('Input  Volume: Tensor of size [2,1,256,256,256]')
    print('Target Volume: Tensor of size [2,1,256,256,256]')
    input_volume  = torch.nn.Parameter(torch.zeros((2,1,256,256,256), device=f'cuda:{PYTORCH_GPU_ID}', requires_grad=True))
    target_volume = torch.ones((2,1,256,256,256), device=f'cuda:{PYTORCH_GPU_ID}', requires_grad=True)
    
    ### Instanciate the geometry
    geo = get_default_3Dgeometry()

    ### Instanciate the differentiable modules
    forward_operator = A(
        geo, np.linspace(0, np.pi, 200), TIGRE_GPU_ID
    )
    backward_operator = At(
        geo, np.linspace(0, np.pi, 200), TIGRE_GPU_ID
    )

    ### Instanciate a simple convolutional layer, the loss and the optimiser
    conv = torch.nn.Conv3d(1,1,1).to(device=torch.device(f'cuda:{PYTORCH_GPU_ID}'))
    mse = torch.nn.MSELoss()
    optimiser = torch.optim.Adam(
        [input_volume] + list(conv.parameters()), 1e-4
    )
    
    ### Compute sinograms
    input_sinogram  = forward_operator(input_volume)
    target_sinogram = forward_operator(target_volume)

    print(f'Testing the forward operator Ax for {n_iterations} iterations')
    for i in range(n_iterations):
        optimiser.zero_grad()
        approximation = conv(forward_operator(input_volume))
        loss = mse(approximation, target_sinogram)
        print(f'MSE Loss: {loss.item():.3f}')
        loss.backward()
        optimiser.step()
    print('Forward operator tests done!')

    print(f'Testing the backward operator Atb for {n_iterations} iterations')
    for i in range(n_iterations):
        optimiser.zero_grad()
        approximation = conv(backward_operator(input_sinogram))
        loss = mse(approximation, target_volume)
        print(f'MSE Loss: {loss.item():.5f}')
        loss.backward()
        optimiser.step()
    print('Backward operator tests done!')


from tigre.utilities.pytorch_bindings import create_pytorch_operator

if __name__ == '__main__':    
    ### Get GPU id
    listGpuNames = gpu.getGpuNames()
    if len(listGpuNames) == 0:
        print("Error: No gpu found")
        raise RuntimeError('There needs to be a GPU')
    else:
        for id in range(len(listGpuNames)):
            print("{}: {}".format(id, listGpuNames[id]))

    TIGRE_GPU_ID   = gpu.getGpuIds(listGpuNames[3])
    TIGRE_GPU_ID = TIGRE_GPU_ID[3]
    PYTORCH_GPU_ID = TIGRE_GPU_ID.devices[0]
    
    print(f'Using GPU {TIGRE_GPU_ID} for TIGRE and GPU {PYTORCH_GPU_ID} for PyTorch')

    geo=get_default_2Dgeometry()
    angles = np.linspace(0, np.pi, 200)
    op, opt = create_pytorch_operator(geo, angles, TIGRE_GPU_ID)
    input_volume = torch.zeros([2,geo.nVoxel[0], geo.nVoxel[1],geo.nVoxel[2]], device=f'cuda:{PYTORCH_GPU_ID}', requires_grad=True)
    sinogram = op(input_volume)
    print(f'Input volume: {input_volume.shape}')
    print(f'Sinogram: {sinogram.shape}')
    print(f'Output volume: {opt(sinogram).shape}')
    #test_2D_operators()
    #test_3D_operators()

    

    