from typing import List
from tigre.utilities.geometry import Geometry
import torch 
import tigre 
import numpy as np

import itertools


def _gpuids(gpuids):
    """None, a GpuIds, or a list of device numbers (int or str) -> None or GpuIds.

    None means: a CUDA tensor runs on its own GPU, data on the host on all GPUs.
    """
    if gpuids is None or hasattr(gpuids, "devices"):
        return gpuids
    from tigre.utilities.gpu import GpuIds
    return GpuIds._from_devices(gpuids)


def _call(op, x, geo, **kwargs):
    """op(x, geo, geo.angles, **kwargs), retrying once if the GPU is out of memory.

    PyTorch keeps the GPU memory it freed in its cache, for its own reuse; TIGRE
    allocates with CUDA directly and cannot use it. If TIGRE runs out of memory,
    give that cache back to CUDA and try again.
    """
    try:
        return op(x, geo, geo.angles, **kwargs)
    except RuntimeError as err:
        if "out of memory" not in str(err) or not torch.cuda.is_initialized():
            raise
    # Outside the except block: by now the failed attempt's exception and its
    # traceback are gone, and with them any GPU tensor the attempt allocated,
    # which would otherwise keep its block of the cache from being freed.
    torch.cuda.empty_cache()
    return op(x, geo, geo.angles, **kwargs)


def _apply(op, x, geo, output_dimension, gpuids, **kwargs):
    """Apply op (tigre.Ax or tigre.Atb) to each 3D array in x, a tensor of
    shape (..., *3D shape). A CUDA tensor stays on the GPU (see tigre.Ax);
    a tensor on the CPU goes through NumPy."""
    gpuids = _gpuids(gpuids)
    x = x.detach()
    extra_dims = list(x.size()[:-3])
    if x.is_cuda and not extra_dims:
        return _call(op, x, geo, gpuids=gpuids, **kwargs)
    output = x.new_empty(extra_dims + list(output_dimension), dtype=torch.float32)
    if not x.is_cuda:
        x = x.numpy()
    for subspace in itertools.product(*[range(dim_size) for dim_size in extra_dims]):
        result = _call(op, x[subspace], geo, gpuids=gpuids, **kwargs)
        output[subspace] = result if torch.is_tensor(result) else torch.from_numpy(result)
    return output

class AXFunction(torch.autograd.Function):
    @staticmethod
    def forward(
        volume_input:torch.Tensor, 
        geo : Geometry,
        gpuids : List[str],
        volume_dimension : List,
        sinogram_dimension : List,
        backprojection_type : str = "matched",
        ) -> torch.Tensor:
        """AXFunction is an autograd function that performs a forward operation
        and preserves the gradients.

        Args:
            volume_input (torch.Tensor): input volume tensor
            geo (Geometry): geometry
            gpuids: GpuIds, list of device numbers, or None (see create_pytorch_operator)
            volume_dimension (List): integers describing the volume dimension
            sinogram_dimension (List): integers describing the volume dimension
            backprojection_type (str): weights of the backprojector, "matched"
                (default, closest to the transpose of Ax) or "FDK"

        Returns:
            torch.Tensor: sinogram
        """
        return _apply(tigre.Ax, volume_input, geo, sinogram_dimension, gpuids)
    
    @staticmethod
    def setup_context(ctx, inputs, output):
        input, geo, gpuids, volume_dimension, sinogram_dimension, backprojection_type = inputs
        ctx.geo = geo
        ctx.volume_dimension = volume_dimension
        ctx.backprojection_type = backprojection_type
        ctx.gpuids = gpuids
        
    @staticmethod
    def backward(ctx, grad_output:torch.Tensor): #type:ignore
        volume_output = _apply(tigre.Atb, grad_output, ctx.geo, ctx.volume_dimension, ctx.gpuids,
                               backprojection_type=ctx.backprojection_type)
        return volume_output, None, None, None, None, None
    
class A(torch.nn.Module):
    """Forward projection as a PyTorch module, with gradients.

    The gradient is computed with TIGRE's backprojector; ``backprojection_type``
    chooses its weights ("matched" by default, the closest to the transpose of
    the forward projector).
    """
    def __init__(self, geo:Geometry, angles, gpuids:List[str], backprojection_type:str="matched"):
        super(A, self).__init__()
        self.backprojection_type = backprojection_type
        geo.angles = angles
        self.geo = geo
        self.gpuids = gpuids
        self.volume_dimension   = geo.nVoxel.tolist()
        self.sinogram_dimension = [len(geo.angles)] + geo.nDetector.tolist()

    def forward(self, x:torch.Tensor):
        return AXFunction.apply(x, self.geo, self.gpuids, self.volume_dimension,
                                self.sinogram_dimension, self.backprojection_type)

class ATBFunction(torch.autograd.Function):
    @staticmethod
    def forward(
        sinogram_input:torch.Tensor, 
        geo : Geometry,
        gpuids : List[str],
        volume_dimension : List,
        sinogram_dimension : List,
        backprojection_type : str = "matched",
        ) -> torch.Tensor:
        """AXFunction is an autograd function that performs a forward operation
        and preserves the gradients.

        Args:
            sinogram_input (torch.Tensor): input sinogram tensor
            geo (Geometry): geometry
            gpuids: GpuIds, list of device numbers, or None (see create_pytorch_operator)
            volume_dimension (List): integers describing the volume dimension
            sinogram_dimension (List): integers describing the volume dimension
            backprojection_type (str): weights of the backprojector, "matched"
                (default, closest to the transpose of Ax) or "FDK"

        Returns:
            torch.Tensor: volume
        """
        return _apply(tigre.Atb, sinogram_input, geo, volume_dimension, gpuids,
                      backprojection_type=backprojection_type)
    
    @staticmethod
    def setup_context(ctx, inputs, output):
        input, geo, gpuids, volume_dimension, sinogram_dimension, backprojection_type = inputs
        ctx.geo = geo
        ctx.sinogram_dimension = sinogram_dimension
        ctx.gpuids = gpuids
        
    @staticmethod
    def backward(ctx, grad_output:torch.Tensor): #type:ignore
        sinogram_output = _apply(tigre.Ax, grad_output, ctx.geo, ctx.sinogram_dimension, ctx.gpuids)
        return sinogram_output, None, None, None, None, None
    
class At(torch.nn.Module):
    """Backprojection as a PyTorch module, with gradients.

    ``backprojection_type`` chooses the backprojector's weights: "matched"
    (default), the closest to the transpose of the forward projector, or
    "FDK", the weights FDK uses.
    """
    def __init__(self, geo:Geometry, angles:np.array, gpuids:List[str], backprojection_type:str="matched"):
        super(At, self).__init__()
        self.backprojection_type = backprojection_type
        geo.angles = angles
        self.geo = geo
        self.gpuids = gpuids
        self.volume_dimension   = geo.nVoxel.tolist()
        self.sinogram_dimension = [len(geo.angles)] + geo.nDetector.tolist()

    def forward(self, x:torch.Tensor):
        return ATBFunction.apply(x, self.geo, self.gpuids, self.volume_dimension,
                                 self.sinogram_dimension, self.backprojection_type)
    




def create_pytorch_operator(geo, angles, gpuids=None, backprojection_type="matched"):
    """Return the forward projector and backprojector as PyTorch modules.

    ``backprojection_type`` sets the weights of the backprojector, used by
    ``At`` and by the gradient of ``A``: "matched" (default, the
    closest to the transpose of the forward projector) or "FDK" (what earlier
    versions used).

    CUDA tensors stay on their GPU: TIGRE reads and writes them there, without
    copies through host memory. ``gpuids`` (a GpuIds or a list of device
    numbers) is only needed to choose GPUs for tensors on the CPU, or to send
    CUDA tensors through the host to several GPUs. Default: a CUDA tensor's
    own GPU, all GPUs for a tensor on the CPU.
    """
    return (A(geo, angles, gpuids, backprojection_type),
            At(geo, angles, gpuids, backprojection_type))



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