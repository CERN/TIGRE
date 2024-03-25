#% DEMO 24: Training a 2D LPD using TIGRE as a backend
#
#
#
#  We interface TIGRE with PyTorch to train a Learned Primal-Dual algorithm
#  to reconstruct MNIST digits
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
# Coded by:           Ander Biguri
# --------------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
plt.style.use('grayscale')

import torch
import torch.nn as nn
from torchvision import datasets, transforms

import tigre
from tigre.utilities.pytorch import Operators
from tigre.utilities import gpu, power_method

def get_default_geometry():
    geo = tigre.geometry()
    #%% Geometry
    geo.DSD = 1536  # Distance Source Detector      (mm) #type:ignore
    geo.DSO = 1000  # Distance Source Origin        (mm) #type:ignore
    # Image parameters
    geo.nVoxel = np.array([1, 28, 28])  # number of voxels              (vx)
    geo.sVoxel = np.array([1, 28, 28])  # total size of the image       (mm)
    geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
    # Detector parameters
    geo.nDetector = np.array([1, 256])  # number of pixels              (px)
    geo.dDetector = np.array([geo.dVoxel[0], 0.8])  # size of each pixel            (mm)
    geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm) #type:ignore
    # Offsets
    geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
    geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
    # MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!
    geo.mode = "parallel"
    return geo

def sinogram_add_noise(
    proj, I0=1000, sigma=5, cross_talk=0.05, flat_field=None, dark_field=None
):
    """
    Adds realistic noise to sinograms.
    - Poisson noise, with I0 counts in a scanner with no sample (bigger value==less noise)
    - Gaussian noise of zero mean and sigma std
    - Detector crosstalk in % of the signal of adjacent pixels.
    - Add a flat_field to add even more realistic noise (computed at non-corrected flat fields)
    """
    dev = torch.cuda.current_device()
    if torch.is_tensor(proj):
        istorch = True
        dev = proj.get_device()
        if dev == -1:
            dev = torch.device("cpu")
    elif isinstance(proj, np.ndarray):
        # all good
        istorch = False
        proj = torch.from_numpy(proj).cuda(dev)
    else:
        raise ValueError("numpy or torch tensor expected")
    if dark_field is None:
        dark_field = torch.zeros(proj.shape, device=dev)
    if flat_field is None:
        flat_field = torch.ones(proj.shape, device=dev)
    max_val = torch.amax(
        proj
    )  # alternatively the highest power of 2 close to this value, but lets leave it as is.

    Im = I0 * torch.exp(-proj / max_val)

    # Uncorrect the flat fields
    Im = Im * (flat_field - dark_field) + dark_field

    # Add Poisson noise
    Im = torch.poisson(Im)

    # Detector cross talk

    kernel = torch.tensor(
        [[0.0, 0.0, 0.0], [cross_talk, 1, cross_talk], [0.0, 0.0, 0.0]]
    ).view(1, 1, 3, 3).repeat(1, 1, 1, 1) / (1 + 2 * cross_talk)

    conv = torch.nn.Conv2d(1, 1, 3, bias=False, padding="same")
    with torch.no_grad():
        conv.weight = torch.nn.Parameter(kernel)
    conv = conv.to(dev)

    Im = conv(Im)
    # Electronic noise:
    Im = Im + sigma * torch.randn(Im.shape, device=dev)

    Im[Im <= 0] = 1e-6
    # Correct flat fields
    Im = (Im - dark_field) / (flat_field - dark_field)
    proj = -torch.log(Im / I0) * max_val
    proj[proj < 0] = 0
    if istorch:
        return proj
    else:
        return proj.cpu().detach().numpy()

def weights_init(module: nn.Module):
    if isinstance(module, nn.Conv2d):
        shape = module.weight.shape
        lim = (6 / (shape[0] + shape[1]) / shape[2] / shape[3]) ** 0.5
        module.weight.data.uniform_(-lim, lim)
        module.bias.data.fill_(0)  # type:ignore

class CnnBlock(nn.Module):
    def __init__(
        self,
        input_channels  = 1,
        output_channels = 1,
        n_filters = 32
    ) -> None:
        super(CnnBlock, self).__init__()

        self.block = nn.Sequential(
            nn.Conv2d(input_channels, n_filters, 3, padding=1),
            nn.PReLU(num_parameters=n_filters, init=0.0),
            nn.Conv2d(n_filters, n_filters, 3, padding=1),
            nn.PReLU(num_parameters=n_filters, init=0.0),
            nn.Conv2d(n_filters, output_channels, 3, padding=1),
        )

        weights_init(self.block)

    def forward(self, input_tensor: torch.Tensor) -> torch.Tensor:
        return self.block(input_tensor)

class Iteration(nn.Module):
    def __init__(
        self,
        A,
        A_T,
        device: torch.device
    ):
        super(Iteration, self).__init__()

        self.A = A
        self.A_T = A_T

        self.operator_norm:float = power_method.svd_power_method(
            arr = np.float32(np.random.random(self.A.geo.nVoxel)),
            geo = self.A.geo,
            angles = self.A.angles
        )

        self.primal_block = CnnBlock(2,1)
        self.primal_block = self.primal_block.to(device)

        self.dual_block = CnnBlock(3,1)
        self.dual_block = self.dual_block.to(device)

    def dual_operation(
        self, primal: torch.Tensor, dual: torch.Tensor, input_sinogram: torch.Tensor
    ) -> torch.Tensor:
        return dual + self.dual_block(
            torch.cat(
                [
                    dual,
                    self.A(primal[:, 0:1, ...]) / self.operator_norm,
                    input_sinogram,
                ],
                dim=1,
            )
        )

    def primal_operation(
        self, primal: torch.Tensor, dual: torch.Tensor
    ) -> torch.Tensor:
        return primal + self.primal_block(
                torch.cat(
                    [
                        primal,
                        self.A_T(dual[:, 0:1, ...]) / self.operator_norm
                    ],
                    dim=1
                )
            )

    def forward(
        self, primal: torch.Tensor, dual: torch.Tensor, input_sinogram: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        # dual block
        dual = self.dual_operation(primal, dual, input_sinogram)
        # primal block
        return self.primal_operation(primal, dual), dual

class LearnedPrimalDual(nn.Module):
    def __init__(
        self,
        A,
        A_T,
        device: torch.device,
        n_iterations = 2
    ):
        super(LearnedPrimalDual, self).__init__()
        self.domain_shape = list(A.geo.nVoxel[1:])
        self.range_shape  = [len(A.angles), A.geo.nDetector[-1]]

        self.device = device
        self.n_iterations = n_iterations
        self.iteration_modules = torch.nn.ModuleDict({
                f"iteration_{i}": Iteration(
                    A,
                    A_T,
                    device
                ) for i in range(self.n_iterations)
            })

    def forward(self, input_sinogram: torch.Tensor,) -> torch.Tensor:

        primal = torch.zeros(
            [input_sinogram.size()[0], 1] + self.domain_shape
            ).to(self.device)

        dual = torch.zeros(
            [input_sinogram.size()[0], 1] + self.range_shape
            ).to(self.device)

        for i in range(self.n_iterations):
            primal, dual = self.iteration_modules[f"iteration_{i}"](
                primal, dual, input_sinogram
            )

        return primal[:, 0:1]

if __name__ == '__main__':

    #%% Define angles of projection and load phantom image
    # define projection angles (in radians)
    angles = np.linspace(0, 2 * np.pi, 256)
     #%% Set GPU ids (we assume only one GPU is available)
    gpu_index = 0
    gpuids = gpu.getGpuIds()
    gpuids.devices = [gpu_index]

    geo = get_default_geometry()

    device = torch.device(f'cuda:{gpu_index}')

    #%% Dataloader boilerplate
    batch_size = 8

    transform=transforms.Compose([
        transforms.ToTensor(),
        transforms.Normalize((0.1307,), (0.3081,))
        ])

    dataset1 = datasets.MNIST('../tigre/data', train=True, download=True,
                       transform=transform)
    train_loader = torch.utils.data.DataLoader(dataset1, batch_size)

    #%% Define torch bindings
    # Perform a L2-loss minimisation of the data discrepancy term (image space)
    loss = torch.nn.MSELoss()

    # Define the operator
    A   = Operators.ForwardOperatorModule(geo, angles, gpuids)
    A_T = Operators.BackwardOperatorModule(geo, angles, gpuids)

    model = LearnedPrimalDual(
        A,
        A_T,
        device,
        n_iterations = 5
    )

    # Define the optimiser
    optimiser = torch.optim.Adam(model.parameters(), lr=1e-4)

    n_steps = 10

    for i in range(n_steps):
        for batch_idx, (data, target) in enumerate(train_loader):
            data = data.to(device)

            with torch.no_grad():
                sinogram = sinogram_add_noise(A(data))

            optimiser.zero_grad()
            rec = model(sinogram)
            mse_loss = loss(rec, data)
            mse_loss.backward()
            optimiser.step()
            print(batch_idx, mse_loss.item())
            if batch_idx % 50 == 0:
                plt.matshow(rec[0,0].detach().cpu())
                plt.savefig('test_lpd.jpg')
                plt.clf()
                plt.close()
