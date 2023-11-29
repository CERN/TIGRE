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
        THREE_D = True
        if len(x.shape) == 4:
            THREE_D = False
        for batch_index in range(x.shape[0]):
            if THREE_D:
                result.append(
                    tigre.Ax(
                        x[batch_index,0],
                        geo,
                        angles,
                        gpuids = gpuids)
                    )
            else:
                result.append(
                    tigre.Ax(
                        x[batch_index],
                        geo,
                        angles,
                        gpuids = gpuids).transpose((1,0,2))
                    )
        if THREE_D:
            return torch.tensor(np.stack(result), requires_grad=True).unsqueeze(1).to(device)
        else:
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