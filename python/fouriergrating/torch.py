"""PyTorch adapter (Plan section "M19: Schnittstellen"): efficiencies as a differentiable function of
numbers of the configuration. Optional: importing this module requires ``torch``.

    import torch
    from fouriergrating.torch import efficiencies

    sim = fg.Simulation(cfg)
    params = ["/layers/0/thickness_um", "/layers/0/breakpoints_um/1"]
    x = torch.tensor([0.3, 0.25], dtype=torch.float64, requires_grad=True)
    e = efficiencies(sim, params, x, [("R", -1, 0)])
    (-e[0]).backward()          # x.grad = -d e_-1 / d x  (discrete adjoint, one VJP)

Forward: the values are written with ``Simulation.set_parameter`` and the efficiencies of ``orders``
(``(side, n, m)``) are returned (0 for evanescent orders); backward: ``Simulation.vjp`` at the same values.
Complex parameters (permittivity) are passed as their real and imaginary parts (``.../n/0``, ``.../n/1``),
which matches the PyTorch convention for real leaf tensors.
"""

from typing import Sequence, Tuple

import torch

from . import Simulation

__all__ = ["efficiencies"]


class _Efficiencies(torch.autograd.Function):
    @staticmethod
    def forward(ctx, values, sim, params, orders):  # type: ignore[override]
        for p, v in zip(params, values.detach().cpu().tolist()):
            sim.set_parameter(p, float(v))
        r = sim.solve()
        out = torch.tensor([r.efficiency(s, n, m) for (s, n, m) in orders], dtype=values.dtype, device=values.device)
        ctx.sim, ctx.params, ctx.orders = sim, list(params), list(orders)
        ctx.values = values.detach().cpu().tolist()
        return out

    @staticmethod
    def backward(ctx, grad_out):  # type: ignore[override]
        for p, v in zip(ctx.params, ctx.values):  # the state of the forward pass
            ctx.sim.set_parameter(p, float(v))
        cot = [(s, n, m, float(g)) for (s, n, m), g in zip(ctx.orders, grad_out.detach().cpu().tolist())]
        g = ctx.sim.vjp(cot, ctx.params)
        return torch.as_tensor(g, dtype=grad_out.dtype, device=grad_out.device), None, None, None


def efficiencies(sim: Simulation, params: Sequence[str], values: "torch.Tensor",
                 orders: Sequence[Tuple[str, int, int]]) -> "torch.Tensor":
    """Efficiencies of ``orders`` for the parameter ``values`` (1D tensor, one entry per pointer)."""
    return _Efficiencies.apply(values, sim, list(params), [(str(s).upper(), int(n), int(m)) for s, n, m in orders])
