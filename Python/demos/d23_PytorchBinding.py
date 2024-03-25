import numpy as np

import tigre
from tigre.utilities.pytorch import Operator
from tigre.utilities import sl3d, gpu

import torch

import matplotlib.pyplot as plt


if __name__ == '__main__':

    #%% Define default geometry
    geo = tigre.geometry_default(high_resolution=False)

    #%% Define angles of projection and load phantom image
    # define projection angles (in radians)
    angles = np.linspace(0, 2 * np.pi, 50)
    # load phantom image
    phantom_type = (
        "yu-ye-wang"  # Default of Python TIGRE Shepp-Logan phantom. Improved visual perception
    )
    shepp_phantom:np.ndarray = sl3d.shepp_logan_3d(
        geo.nVoxel, phantom_type=phantom_type
    )
    print(f'Phantom shape: {shepp_phantom.shape}')

    #%% Set GPU ids (we assume only one GPU is available)
    gpu_index = 2
    gpuids = gpu.getGpuIds()
    gpuids.devices = [gpu_index]

    #%% Define torch bindings
    # Perform a L2-loss minimisation of the data discrepancy term
    loss = torch.nn.MSELoss()
    x = torch.zeros(
        size = (1,1)+shepp_phantom.shape,
        device=f'cuda:{gpu_index}',
        requires_grad=True
        )
    x_target = torch.from_numpy(shepp_phantom).to(f'cuda:{gpu_index}').unsqueeze(0).unsqueeze(0)
    assert x.size() == x_target.size()

    # Define the operator
    operator = Operator.OperatorModule(geo, angles, gpuids)

    # Define the optimiser
    optimiser = torch.optim.SGD([x], lr=1)

    n_steps = 101

    for i in range(n_steps):
        optimiser.zero_grad()
        y = operator(x)
        mse_loss = loss(y, operator(x_target))
        mse_loss.backward()
        optimiser.step()

    plt.style.use('grayscale')
    plt.matshow(x_target[0,0,32].detach().cpu())
    plt.savefig(f'target.jpg')
    plt.clf()
    plt.matshow(x[0,0,32].detach().cpu())
    plt.savefig(f'reconstruction.jpg')
    plt.clf()


