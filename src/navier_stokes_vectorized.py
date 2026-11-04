"""Vectorized explicit momentum step for the Navier-Stokes solver.

These functions reproduce the per-cell finite-volume stencils of
``navier_stokes.py`` (``navier_stokes_x`` / ``navier_stokes_y``) using NumPy
array slicing, so the whole interior is updated in one shot. The arithmetic is
identical to the scalar version, only the iteration is removed.
"""

import numpy as np


def navier_stokes_x(u, v, dx, dy, dt, Re, out=None):
    """x-momentum update on the staggered grid -> returns ``u_tp``.

    Interior nodes are i = 1..Nxc-1, j = 1..Nyc (matching the scalar loops
    ``for i in range(1, Nxc)`` / ``for j in range(1, Nyc+1)``).

    Pass ``out`` to write into a preallocated buffer. Cells this stencil does
    not touch are left as they are, so ``out`` must start at zero.
    """

    u_tp = np.zeros_like(u) if out is None else out
    nx, ny = u.shape

    i0, i1 = 1, nx - 1   # i = 1 .. nx-2  (= Nxc-1)
    j0, j1 = 1, ny - 1   # j = 1 .. ny-2  (= Nyc)

    # d(uu)/dx
    uu_f = (0.5 * (u[i0+1:i1+1, j0:j1] + u[i0:i1, j0:j1]))**2.0
    uu_b = (0.5 * (u[i0-1:i1-1, j0:j1] + u[i0:i1, j0:j1]))**2.0
    Duu_x = (uu_f - uu_b) / dx

    # d(uv)/dy
    u_sum_f = u[i0:i1, j0:j1] + u[i0:i1, j0+1:j1+1]
    u_sum_b = u[i0:i1, j0:j1] + u[i0:i1, j0-1:j1-1]
    v_sum_f = v[i0:i1, j0:j1] + v[i0+1:i1+1, j0:j1]
    v_sum_b = v[i0:i1, j0-1:j1-1] + v[i0+1:i1+1, j0-1:j1-1]
    uv_x_f = 0.25 * u_sum_f * v_sum_f
    uv_x_b = 0.25 * u_sum_b * v_sum_b
    Duv_y = (uv_x_f - uv_x_b) / dy

    # Laplacian of u
    D2u_x2 = (u[i0+1:i1+1, j0:j1] - 2.0*u[i0:i1, j0:j1] + u[i0-1:i1-1, j0:j1]) / dx**2.0
    D2u_y2 = (u[i0:i1, j0+1:j1+1] - 2.0*u[i0:i1, j0:j1] + u[i0:i1, j0-1:j1-1]) / dy**2.0

    advection_u = Duu_x + Duv_y
    diffusion_u = D2u_x2 + D2u_y2

    u_tp[i0:i1, j0:j1] = u[i0:i1, j0:j1] + dt * (-advection_u + (1.0/Re) * diffusion_u)

    return u_tp


def navier_stokes_y(u, v, dx, dy, dt, Re, out=None):
    """y-momentum update on the staggered grid -> returns ``v_tp``.

    Interior nodes are i = 1..Nxc, j = 1..Nyc-1 (matching the scalar loops
    ``for i in range(1, Nxc+1)`` / ``for j in range(1, Nyc)``).

    Pass ``out`` to write into a preallocated buffer. Cells this stencil does
    not touch are left as they are, so ``out`` must start at zero.
    """

    v_tp = np.zeros_like(v) if out is None else out
    nx, ny = v.shape

    i0, i1 = 1, nx - 1   # i = 1 .. nx-2  (= Nxc)
    j0, j1 = 1, ny - 1   # j = 1 .. ny-2  (= Nyc-1)

    # d(vv)/dy
    vv_f = (0.5 * (v[i0:i1, j0+1:j1+1] + v[i0:i1, j0:j1]))**2.0
    vv_b = (0.5 * (v[i0:i1, j0:j1] + v[i0:i1, j0-1:j1-1]))**2.0
    Dvv_y = (vv_f - vv_b) / dy

    # d(uv)/dx
    u_sum_f = u[i0:i1, j0+1:j1+1] + u[i0:i1, j0:j1]
    u_sum_b = u[i0-1:i1-1, j0+1:j1+1] + u[i0-1:i1-1, j0:j1]
    v_sum_f = v[i0:i1, j0:j1] + v[i0+1:i1+1, j0:j1]
    v_sum_b = v[i0:i1, j0:j1] + v[i0-1:i1-1, j0:j1]
    uv_y_f = 0.25 * u_sum_f * v_sum_f
    uv_y_b = 0.25 * u_sum_b * v_sum_b
    Duv_x = (uv_y_f - uv_y_b) / dx

    # Laplacian of v
    D2v_x2 = (v[i0+1:i1+1, j0:j1] - 2.0*v[i0:i1, j0:j1] + v[i0-1:i1-1, j0:j1]) / dx**2.0
    D2v_y2 = (v[i0:i1, j0+1:j1+1] - 2.0*v[i0:i1, j0:j1] + v[i0:i1, j0-1:j1-1]) / dy**2.0

    advection_v = Duv_x + Dvv_y
    diffusion_v = D2v_x2 + D2v_y2

    v_tp[i0:i1, j0:j1] = v[i0:i1, j0:j1] + dt * (-advection_v + (1.0/Re) * diffusion_v)

    return v_tp
