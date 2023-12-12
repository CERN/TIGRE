from typing import List

import tigre
import torch
import numpy as np

### Is there a TIGRE geo class?

class OperatorFunction(torch.autograd.Function):
    @staticmethod
    def forward(x:torch.Tensor, geo, angles:np.ndarray, gpuids:List) -> torch.Tensor:
        
        # We will allow 2D and 3D (i.e. 4D and 5D tensors)
        # The code will work for 5D tensors, so if its 4D, lets make it 5D
        if geo.nVoxel[0]==1: # i.e. if len(x.shape)==4
            if len(x.shape)!=4:
                raise ValueError("Dimensions of tensor don't match expected dimensions for a 2D CT case")
            x.unsqueeze(2) 
        else:
            if len(x.shape)!=5:
                raise ValueError("Dimensions of tensor don't match expected dimensions for a 3D CT case")
            
        # For now we only support channel-by-channel. 
        if x.shape(1)!=1:
            raise NotImplementedError("TIGRE torch operator only accepts 1 channel (for now). Contact the devs if you have a reason to support more")

        device = x.device
        x = x.detach().cpu().numpy()

        result = np.zeros((x.shape[0],x.shape[1],len(angles)))
        
        for batch_index in range(x.shape[0]):
            result.append(
                tigre.Ax(
                    x[batch_index,0],
                    geo,
                    angles,
                    gpuids = gpuids)
                )
            
        # The return will be different depending if its 2D or 3D tomography too. 
        if geo.nVoxel[0]==1: # i.e. if len(x.shape)==4
            return torch.tensor(result, requires_grad=True).squeeze(2).to(device)
        else:
            return torch.tensor(result, requires_grad=True).to(device)

    @staticmethod
    def setup_context(ctx, inputs, output) -> None:
        input, geo, angles, gpuids = inputs
        ctx.geo = geo
        ctx.angles = angles
        ctx.gpuids = gpuids

    @staticmethod
    def backward(ctx, grad_output):
        device = grad_output.device
        THREE_D = True
        if len(grad_output.shape) == 4:
            THREE_D = False
        if THREE_D:
            print(grad_output.size())
            grad_output = grad_output.detach().cpu().numpy().transpose((0,1,3,2,4))
        else:
            grad_output = grad_output.detach().cpu().numpy().transpose((0,2,1,3))
        result = []

        for batch_index in range(grad_output.shape[0]):
            if THREE_D:
                result.append(
                    tigre.Atb(
                        grad_output[batch_index][0],
                        ctx.geo,
                        ctx.angles,
                        gpuids = ctx.gpuids)
                    )
            else:
                result.append(
                    tigre.Atb(
                        grad_output[batch_index],
                        ctx.geo,
                        ctx.angles,
                        gpuids = ctx.gpuids)
                    )

        if THREE_D:
            return torch.tensor(np.stack(result), requires_grad=True).unsqueeze(1).to(device), None, None, None
        else:
            return torch.tensor(np.stack(result), requires_grad=True).to(device), None, None, None


class OperatorModule(torch.nn.Module):
    def __init__(self, geo, angles, gpuids: List):
        super(OperatorModule, self).__init__()
        self.geo = geo
        self.angles = angles
        self.gpuids = gpuids

    def forward(self, x:torch.Tensor):
        return OperatorFunction.apply(x, self.geo, self.angles, self.gpuids)