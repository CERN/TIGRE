import tigre

import torch
import numpy as np

class AXFunction(torch.autograd.Function):
    @staticmethod
    def forward(
        input_tensor:torch.Tensor, 
        geo, 
        angles, 
        gpuids
        ):
        input_dimension = input_tensor.size()
        device = input_tensor.device
        input_tensor = input_tensor.detach().cpu().numpy()
        assert input_dimension[1] == 1, f'The input tensor has a its first dimension equal to {input_dimension[1]}. Currently, only tensors with channel dimension = 1 are supported.'
        if len(input_dimension) == 5:
            input_tensor = input_tensor.squeeze(1)
        result = []
        for batch_index in range(input_tensor.shape[0]):
            ### This is a hack.
            ### TIGRE expects an array of shape [NVoxDepth, NVoxHeight, NVoxWidth]
            ### Pytorch provides an array of shape 
            #   -> [BatchSize, NChannels = 1, NVoxHeight, NVoxWidth] for 2D Tomography
            #   -> [BatchSize, NChannels = 1, NVoxDepth, NVoxHeight, NVoxWidth] for 3D Tomography
            # How we handle 2D case: we assume the channel dimension (=1) for PyTorch is the NVoxDepth (=1) dimension for TIGRE
            # How we handle 3D case: we squeeze the input tensor to remove the channels dimension 
            ax:np.ndarray = tigre.Ax(
                    input_tensor[batch_index], 
                    geo,
                    angles,
                    gpuids = gpuids)
            # Ax returns an array of size [N_angles, DetX, DetY]
            # For 2D, DetX = 1, and we want to output a tensor of size [BatchSize, NChannels = 1, N_angles, DetY]
            # We then use the fact that DetX = 1 and transpose the output array
            if len(input_dimension) == 4:
                ax = ax.transpose((1,0,2))
            result.append(
                ax
                )
        result = torch.tensor(np.stack(result), requires_grad=True).to(device)
        if len(input_dimension) == 5:
            result = result.unsqueeze(1)
        return result

    @staticmethod
    def setup_context(ctx, inputs, output):
        input, geo, angles, gpuids = inputs
        ctx.geo = geo
        ctx.angles = angles
        ctx.gpuids = gpuids

    @staticmethod
    def backward(ctx, grad_output:torch.Tensor):
        device = grad_output.device
        input_dimension = grad_output.size()
        grad_output : np.ndarray= grad_output.detach().cpu().numpy()
        if len(input_dimension) == 5:
            grad_output = grad_output.squeeze(1)
        elif len(input_dimension) == 4:
            grad_output = np.expand_dims(grad_output.squeeze(1),axis=2)
            
        result = []
        for batch_index in range(grad_output.shape[0]):
            result[batch_index] = tigre.Atb(
                    grad_output[batch_index],
                    ctx.geo,
                    ctx.angles,
                    gpuids = ctx.gpuids)
        result = torch.tensor(result, requires_grad=True,device=device)

        if len(input_dimension) == 5:
            result = result.unsqueeze(1)
        return result, None, None, None
    
class ATBFunction(torch.autograd.Function):
    @staticmethod
    def forward(input_tensor:torch.Tensor, geo, angles, gpuids):
        device = input_tensor.device
        input_dimension = input_tensor.size()
        input_shape = input_tensor.shape
        input_tensor:np.ndarray = input_tensor.detach().cpu().numpy()
        if len(input_dimension) == 5:
            input_tensor = input_tensor.squeeze(1)
        elif len(input_dimension) == 4:
            input_tensor = np.expand_dims(input_tensor.squeeze(1),axis=2)
            
        result = np.zeros((input_tensor.shape[0],geo.nVoxel[1],geo.nVoxel[2]))
        for batch_index in range(input_tensor.shape[0]):
            result[batch_index]  = tigre.Atb(
                    input_tensor[batch_index],
                    geo,
                    angles,
                    gpuids = gpuids)
        
        result = torch.tensor(result, requires_grad=True,device=device)

        if len(input_dimension) == 5:
            result = result.unsqueeze(1)

        return result

    @staticmethod
    def setup_context(ctx, inputs, output):
        input, geo, angles, gpuids = inputs
        ctx.geo = geo
        ctx.angles = angles
        ctx.gpuids = gpuids

    @staticmethod
    def backward(ctx, grad_output:torch.Tensor):
        input_dimension = grad_output.size()
        device = grad_output.device
        grad_output = grad_output.detach().cpu().numpy()
        assert input_dimension[1] == 1, f'The input tensor has a its first dimension equal to {input_dimension[1]}. Currently, only tensors with channel dimension = 1 are supported.'
        if len(input_dimension) == 5:
            grad_output = grad_output.squeeze(1)
        result = []
        for batch_index in range(grad_output.shape[0]):
            ### This is a hack.
            ### TIGRE expects an array of shape [NVoxDepth, NVoxHeight, NVoxWidth]
            ### Pytorch provides an array of shape 
            #   -> [BatchSize, NChannels = 1, NVoxHeight, NVoxWidth] for 2D Tomography
            #   -> [BatchSize, NChannels = 1, NVoxDepth, NVoxHeight, NVoxWidth] for 3D Tomography
            # How we handle 2D case: we assume the channel dimension (=1) for PyTorch is the NVoxDepth (=1) dimension for TIGRE
            # How we handle 3D case: we squeeze the input tensor to remove the channels dimension 
            ax:np.ndarray = tigre.Ax(
                    grad_output[batch_index], 
                    ctx.geo,
                    ctx.angles,
                    gpuids = ctx.gpuids)
            # Ax returns an array of size [N_angles, DetX, DetY]
            # For 2D, DetX = 1, and we want to output a tensor of size [BatchSize, NChannels = 1, N_angles, DetY]
            # We then use the fact that DetX = 1 and transpose the output array
            if len(input_dimension) == 4:
                ax = ax.transpose((1,0,2))
            result.append(ax)
        result = torch.tensor(np.stack(result), requires_grad=True).to(device)
        if len(input_dimension) == 5:
            result = result.unsqueeze(1)
        return result, None, None, None

class A(torch.nn.Module):
    def __init__(self, geo, angles, gpuids):
        super(A, self).__init__()
        self.geo = geo
        self.angles = angles
        self.gpuids = gpuids

    def forward(self, x):
        return AXFunction.apply(x, self.geo, self.angles, self.gpuids)
    
class At(torch.nn.Module):
    def __init__(self, geo, angles, gpuids):
        super(At, self).__init__()
        self.geo = geo
        self.angles = angles
        self.gpuids = gpuids

    def forward(self, x):
        return ATBFunction.apply(x, self.geo, self.angles, self.gpuids)


def create_pytorch_operator(geo, angles, gpuids):
    return A(geo, angles, gpuids), At(geo, angles, gpuids)

## This may be useful for non-torch stuff, but doesn't work for torch autograd. 
#  I'll leave it here for now.
class Operator:
    def __init__(self, geo, angles, gpuids):
        super(Operator, self).__init__()
        self.geo = geo
        self.angles = angles
        self.gpuids = gpuids
        self.ax = A(self.geo, self.angles, self.gpuids)
        self.atb = At(self.geo, self.angles, self.gpuids)

    def __call__(self, x):
        return self.forward(x)
    def forward(self, x):
        return self.ax(x)
    def T(self,b):
        return self.backward(b)
    def backward(self, b):
        return self.atb(b)