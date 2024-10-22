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
import tigre.algorithms as algs

import numpy as np
import tigre.utilities.gpu as gpu
import matplotlib.pyplot as plt



## Imports:
# You can import torch compatible and autograd-compatible forward and backprojectors as:
from tigre.utilities.pytorch_bindings import A, At
# You can also import the function to create the forward and backprojectors in onle line:
from tigre.utilities.pytorch_bindings import create_pytorch_operator

## Create a 2D geometry
geo = tigre.geometry(mode="fan")
angles = np.linspace(0, np.pi, 200)
# you can now just create operators that are pytorch compatible as:
ax, atb = create_pytorch_operator(geo, angles)
# You may however want to have some level of control on which GPU TIGRE will be using, which you can just define as:
tigre_devices   = gpu.getGpuIds()
tigre_devices = tigre_devices[0]
ax, atb = create_pytorch_operator(geo, angles, tigre_devices)

## Now we can create some torch tensors:
torch_devices = 3
torch.cuda.set_device(torch_devices)
print(f'Using GPU {tigre_devices} for TIGRE and GPU {torch_devices} for PyTorch')

# Note: TIGRE pytorch operators take a batched tensor, i.e. a 4D tensor for 2D tomography and a 5D tensor for 3D tomography.
input_volume = torch.zeros([2,geo.nVoxel[0], geo.nVoxel[1],geo.nVoxel[2]], device=f'cuda:{torch_devices}', requires_grad=True)

sinogram = ax(input_volume)
print(f'Input volume: {input_volume.shape}')
print(f'Sinogram: {sinogram.shape}')
print(f'Output volume: {atb(sinogram).shape}')

# =============================================================================
## Example of use inside a neural network

## We will use Learned Gradient from:
# Solving ill-posed inverse problems using iterative deep neural networks. Adler, Jonas, and Öktem, Ozan. Inverse Problems 2017

# Lets define the learned blocks first. This is just 4 convolutional blocks with ReLUs. 
class LGblock(torch.nn.Module):
    def __init__(self, channels=[8,32,32,6]):
        super(LGblock, self).__init__()
        layers = len(channels) - 1
        layer_list = []
        for ii in range(layers):
            layer_list.append(torch.nn.Conv2d(channels[ii], channels[ii + 1], 3, padding=1))
            # Have ReLUs all the way except the last
            if ii < layers - 1:
                layer_list.append(torch.nn.ReLU())
        self.block = torch.nn.Sequential(*layer_list)
    def forward(self, x):
        return self.block(x)
# Now lets make a LG model
class LearnedGradient(torch.nn.Module):
    def __init__(self, geo, angles, tigre_gpuids=None):
        super(LearnedGradient, self).__init__()
        self.tigre_gpuids = tigre_gpuids
        self.geo = geo
        self.n_iters = 5
        self.ax, self.atb = create_pytorch_operator(geo, angles, self.tigre_gpuids)
        gd_step = 1e-4
        self.M = 6 # last channel size 
        for i in range(self.n_iters):
            self.add_module(f"{i}_conv", LGblock())

        self.step_size = torch.nn.ParameterList(
            [torch.nn.Parameter(torch.ones(1) * gd_step) for i in range(self.n_iters)]
        )
    def forward(self, x):

        B, C, W, H = x.shape
        # do an initial reconstruction
        f = x.new_zeros(B, 1, *self.geo.nVoxel[1:])
        for i in range(B):
            sino = x[i].detach().cpu().numpy()
            aux = algs.fdk(sino,self.geo, angles, gpuids=self.tigre_gpuids)
            aux = torch.clip(torch.from_numpy(aux), min=0).to(x.device)
            f[i] = aux
        # Now we will do the iterations
        tmp = torch.zeros(f.shape).type_as(f)
        s = tmp.clone()
        for i in range(self.M - 1):
            s = torch.cat((s, tmp), dim=1)

        for i in range(self.n_iters):
            conv_module = getattr(self, f"{i}_conv")
            del_L = self.step_size[i] * self.atb(x - self.ax(f))
            output = conv_module(torch.cat([f, s, del_L], dim=1))
            # output last channel
            f = f + output[:, self.M-1 : self.M]
            s = torch.nn.ReLU()(output[:, 0 : self.M])
        return f





# =============================================================================
# We will now train Learned Gradient on MNIST


tigre_devices   = gpu.getGpuIds()
tigre_devices = tigre_devices[3]

## Get dataset:
import torchvision


datapath = '/local/scratch/public/ab2860/data/mnist'
dataset = torchvision.datasets.MNIST(datapath, download=True,transform=torchvision.transforms.PILToTensor())
# lets make it smaller. 
dataset = torch.utils.data.Subset(dataset,  list(range(100)))

# We are going to be reconstructing 28x28 images, so we need to define the geometry accordingly
geo = tigre.geometry(mode="fan",nVoxel=[1,28,28]) # MNIST is 28x28
angles = np.linspace(0, 2*np.pi, 50)

# take one sample
sample_image = dataset[0][0]
sample_image = sample_image.float().numpy()
print(sample_image.shape)
print(geo)
# make sinogram and recon with fdk
sample_sinogram = tigre.Ax(sample_image, geo, angles,gpuids=tigre_devices)
recon = tigre.algorithms.single_pass_algorithms.fdk(sample_sinogram, geo, angles, gpuids=tigre_devices)
# plot
plt.figure()
plt.subplot(131)
plt.imshow(sample_image[0])
plt.title('Input image')
plt.subplot(132)
plt.imshow(sample_sinogram[:,0,:])
plt.title('Sinogram')
plt.subplot(133)
plt.imshow(recon[0])
plt.title('Reconstruction')
plt.savefig('recon_with_fdk.png')

# =============================================================================
## We will now do the training
# lets make a DataLoader
dataloader = torch.utils.data.DataLoader(dataset, batch_size=5, shuffle=True)
# Lets make operators to crate the sinograms:
ax, atb = create_pytorch_operator(geo, angles, tigre_devices)

# Create model
model = LearnedGradient(geo, angles, tigre_devices).to(f'cuda:{torch_devices}')
model.train()
# Create optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
# loss function
loss_fn = torch.nn.MSELoss()

for epoch in range(10):
    for i, (data, _) in enumerate(dataloader):
        # make sinograms for the current batch:
        data = data.float()
        data = data.to(f'cuda:{torch_devices}')
        sinogram = ax(data)
        # zero the gradients
        optimizer.zero_grad()
        # forward pass
        output = model(sinogram)
        # compute loss
        loss = loss_fn(output, data)
        # backward pass
        loss.backward()
        # update weights
        optimizer.step()
    print(f'Loss: {loss.item()}')


# =============================================================================
# Lets see the results
model.eval()
lg_recon = model(torch.from_numpy(sample_sinogram).unsqueeze(0).to(f'cuda:{torch_devices}'))
plt.figure()
plt.subplot(131)
plt.imshow(sample_image[0])
plt.title('Input image')
plt.subplot(132)
plt.imshow(recon[0])
plt.title('FDK')
plt.subplot(133)
plt.imshow(lg_recon[0,0].detach().cpu().numpy())
plt.title('Learned Gradient')
plt.savefig('recon_pytorch.png')
