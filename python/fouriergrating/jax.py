"""JAX adapter (Plan section "M19: Schnittstellen"): efficiencies as a differentiable JAX function of
numbers of the configuration (``jax.custom_vjp`` with ``jax.pure_callback``). Optional: requires ``jax``;
use float64 (``jax.config.update("jax_enable_x64", True)``).

    from fouriergrating.jax import efficiency_function

    f = efficiency_function(fg.Simulation(cfg), ["/layers/0/thickness_um"], [("R", -1, 0)])
    value, grad = jax.value_and_grad(lambda x: -f(x)[0])(jnp.array([0.3]))
"""

from typing import Callable, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

from . import Simulation

__all__ = ["efficiency_function"]


def efficiency_function(sim: Simulation, params: Sequence[str],
                        orders: Sequence[Tuple[str, int, int]]) -> Callable:
    """Returns f(values) -> efficiencies of ``orders`` (``(side, n, m)``) with a custom VJP (discrete adjoint)."""
    params = list(params)
    orders = [(str(s).upper(), int(n), int(m)) for s, n, m in orders]

    def _set(values) -> None:
        for p, v in zip(params, np.asarray(values, dtype=float).tolist()):
            sim.set_parameter(p, float(v))

    def _forward(values):
        _set(values)
        r = sim.solve()
        return np.array([r.efficiency(s, n, m) for (s, n, m) in orders], dtype=np.float64)

    def _vjp(values, g):
        _set(values)
        cot = [(s, n, m, float(w)) for (s, n, m), w in zip(orders, np.asarray(g, dtype=float).tolist())]
        return np.asarray(sim.vjp(cot, params), dtype=np.float64)

    out_shape = jax.ShapeDtypeStruct((len(orders),), jnp.float64)

    @jax.custom_vjp
    def f(values):
        return jax.pure_callback(_forward, out_shape, values)

    def f_fwd(values):
        return f(values), values

    def f_bwd(values, g):
        grad = jax.pure_callback(_vjp, jax.ShapeDtypeStruct(values.shape, jnp.float64), values, g)
        return (grad,)

    f.defvjp(f_fwd, f_bwd)
    return f
