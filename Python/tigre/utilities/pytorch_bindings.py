from typing import List
from tigre.utilities.geometry import Geometry
import torch 
import tigre 
import numpy as np

import itertools

class AXFunction(torch.autograd.Function):
    @staticmethod
    def forward(
        volume_input:torch.Tensor, 
        geo : Geometry,
        gpuids : List[str],
        volume_dimension : List,
        sinogram_dimension : List
        ) -> torch.Tensor:
        """AXFunction is an autograd function that performs a forward operation
        and preserves the gradients.

        Args:
            volume_input (torch.Tensor): input volume tensor
            geo (Geometry): geometry
            gpuids (List[str]): gpus on which to perform the computations
            volume_dimension (List): integers describing the volume dimension
            sinogram_dimension (List): integers describing the volume dimension

        Returns:
            torch.Tensor: sinogram
        """
        device = volume_input.device
        
        ### Might be good to add a check
        extra_dimensions = len(volume_input.size()) - 3
        
        extra_dims = list(volume_input.size()[:extra_dimensions])
        sinogram_output = volume_input.new_empty((extra_dims + sinogram_dimension), dtype=torch.float32) # type:ignore

        volume_input = volume_input.detach().cpu().numpy()  
        for subspace in itertools.product(*[range(dim_size) for dim_size in extra_dims]):
            sinogram_output[subspace] = torch.from_numpy(
                tigre.Ax(volume_input[subspace],geo, geo.angles, gpuids = gpuids)).to(device)
        return sinogram_output
    
    @staticmethod
    def setup_context(ctx, inputs, output):
        input, geo, gpuids, volume_dimension, sinogram_dimension = inputs
        ctx.geo = geo
        ctx.volume_dimension = volume_dimension
        ctx.gpuids = gpuids
        
    @staticmethod
    def backward(ctx, grad_output:torch.Tensor): #type:ignore
        device = grad_output.device
        
        extra_dimensions = len(grad_output.size()) - 3
        
        extra_dims = list(grad_output.size()[:extra_dimensions])
        volume_output = grad_output.new_empty((extra_dims + ctx.volume_dimension), dtype=torch.float32) # type:ignore
        
        grad_output : np.ndarray = grad_output.detach().cpu().numpy()        
        for subspace in itertools.product(*[range(dim_size) for dim_size in extra_dims]):
            volume_output[subspace] = torch.from_numpy(
                tigre.Atb(grad_output[subspace], ctx.geo, ctx.geo.angles, gpuids = ctx.gpuids)).to(device)
        
        return volume_output, None, None, None, None
    
class A(torch.nn.Module):
    def __init__(self, geo:Geometry, angles, gpuids:List[str]):
        super(A, self).__init__()
        geo.angles = angles
        self.geo = geo
        self.gpuids = gpuids
        self.volume_dimension   = geo.nVoxel.tolist()
        self.sinogram_dimension = [len(geo.angles)] + geo.nDetector.tolist()

    def forward(self, x:torch.Tensor):
        return AXFunction.apply(x, self.geo, self.gpuids, self.volume_dimension, self.sinogram_dimension)

class ATBFunction(torch.autograd.Function):
    @staticmethod
    def forward(
        sinogram_input:torch.Tensor, 
        geo : Geometry,
        gpuids : List[str],
        volume_dimension : List,
        sinogram_dimension : List
        ) -> torch.Tensor:
        """AXFunction is an autograd function that performs a forward operation
        and preserves the gradients.

        Args:
            sinogram_input (torch.Tensor): input sinogram tensor
            geo (Geometry): geometry
            gpuids (List[str]): gpus on which to perform the computations
            volume_dimension (List): integers describing the volume dimension
            sinogram_dimension (List): integers describing the volume dimension

        Returns:
            torch.Tensor: volume
        """
        device = sinogram_input.device
        
        ### Might be good to add a check
        extra_dimensions = len(sinogram_input.size()) - 3
        
        extra_dims = list(sinogram_input.size()[:extra_dimensions])
        volume_output = sinogram_input.new_empty((extra_dims + volume_dimension), dtype=torch.float32) # type:ignore

        sinogram_input = sinogram_input.detach().cpu().numpy()  
        for subspace in itertools.product(*[range(dim_size) for dim_size in extra_dims]):
            volume_output[subspace] = torch.from_numpy(
                tigre.Atb(sinogram_input[subspace],geo, geo.angles, gpuids = gpuids)).to(device)
        return volume_output
    
    @staticmethod
    def setup_context(ctx, inputs, output):
        input, geo, gpuids, volume_dimension, sinogram_dimension = inputs
        ctx.geo = geo
        ctx.sinogram_dimension = sinogram_dimension
        ctx.gpuids = gpuids
        
    @staticmethod
    def backward(ctx, grad_output:torch.Tensor): #type:ignore
        device = grad_output.device
        
        extra_dimensions = len(grad_output.size()) - 3
        
        extra_dims = list(grad_output.size()[:extra_dimensions])
        sinogram_output = grad_output.new_empty((extra_dims + ctx.sinogram_dimension), dtype=torch.float32) # type:ignore
        
        grad_output : np.ndarray = grad_output.detach().cpu().numpy()        
        for subspace in itertools.product(*[range(dim_size) for dim_size in extra_dims]):
            sinogram_output[subspace] = torch.from_numpy(
                tigre.Ax(grad_output[subspace], ctx.geo, ctx.geo.angles, gpuids = ctx.gpuids)).to(device)
        
        return sinogram_output, None, None, None, None
    
class At(torch.nn.Module):
    def __init__(self, geo:Geometry, angles:np.array, gpuids:List[str]):
        super(At, self).__init__()
        geo.angles = angles
        self.geo = geo
        self.gpuids = gpuids
        self.volume_dimension   = geo.nVoxel.tolist()
        self.sinogram_dimension = [len(geo.angles)] + geo.nDetector.tolist()

    def forward(self, x:torch.Tensor):
        return ATBFunction.apply(x, self.geo, self.gpuids, self.volume_dimension, self.sinogram_dimension)
    




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