from typing import List

import tigre
import torch
import numpy as np

### Is there a TIGRE geo class?

class OperatorFunction(torch.autograd.Function):
    @staticmethod
    def forward(x:torch.Tensor, geo, angles:np.ndarray, gpuids:List) -> torch.Tensor:
        device = x.device
        x = x.detach().cpu().numpy()
        result = []
        for batch_index in range(x.shape[0]):
            result.append(
                tigre.Ax(
                    x[batch_index],
                    geo,
                    angles,
                    projection_type='FDK',
                    gpuids = gpuids).transpose((1,0,2))
                )
        return torch.tensor(np.stack(result), requires_grad=True).to(device)

    @staticmethod
    def setup_context(ctx, inputs, output) -> None:
        input, geo, angles, gpuids = inputs
        ctx.geo = geo
        ctx.angles = angles
        ctx.gpuids = gpuids

    @staticmethod
    def backward(ctx, grad_output):
        device = grad_output.device
        grad_output = grad_output.detach().cpu().numpy().transpose((0,2,1,3))
        result = []
        for batch_index in range(grad_output.shape[0]):
            result.append(
                tigre.Atb(
                    grad_output[batch_index],
                    ctx.geo,
                    ctx.angles,
                    projection_type='FDK',
                    gpuids = ctx.gpuids)
                )
        return torch.tensor(np.stack(result), requires_grad=True).to(device), None, None, None, None

class OperatorModule(torch.nn.Module):
    def __init__(self, geo, angles, gpuids: List, projection_type:str):
        super(OperatorModule, self).__init__()
        self.geo = geo
        self.angles = angles
        self.gpuids = gpuids
        self.projection_type = projection_type

    def forward(self, x):
        return OperatorFunction.apply(x, self.geo, self.angles, self.gpuids, projection_type=self.projection_type)