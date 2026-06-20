"""Fully vectorized lid-driven cavity Navier-Stokes solver.

Drop-in replacement for ``flow_simulation.simulate_flow`` with an identical
signature and identical return values, but every per-cell Python loop has been
replaced by NumPy array operations:

* velocity / pressure / vorticity boundary conditions -> slice assignments
* explicit momentum update                            -> stencil slicing
* pressure & streamfunction SOR sweeps                -> red-black SOR
* velocity correction and contour fields              -> stencil slicing

The momentum, correction and contour steps are arithmetically identical to the
scalar code; the SOR sweeps use red-black ordering (see ``solver_vectorized``)
and converge to the same flow field.

The explicit time advancement is only conditionally stable. ``simulate_flow``
checks the CFL / diffusive timestep limits up front and warns if ``dt`` is too
large, and it guards against divergence: if the velocity field becomes
non-finite it stops, reports the time at which it blew up, and returns the last
finite state instead of silently producing NaN (blank) charts. Pass
``return_diagnostics=True`` to also get the per-step convergence history.
"""

import warnings

import numpy as np

from navier_stokes_vectorized import navier_stokes_x, navier_stokes_y
from solver_vectorized import poisson, streamlines, checkerboard_masks
from boundary_conditions_vectorized import (
    velocity_bc, pressure_bc, pressure_ghosts, vorticity_bc
)


def _check_stability(dt, dx, dy, Re, Us, Un, Vw, Ve):
    """Warn if ``dt`` violates the explicit-scheme stability limits.

    Returns the convective CFL number and the diffusive timestep limit so they
    can be reported in the diagnostics.
    """

    h = min(dx, dy)
    u_ref = max(abs(Us), abs(Un), abs(Vw), abs(Ve), 1e-30)
    cfl = u_ref * dt / h                 # convective CFL number
    dt_visc = Re * h**2 / 4.0            # diffusive stability limit

    if cfl > 1.0:
        warnings.warn(
            f"Convective CFL = {cfl:.2f} > 1: dt={dt:g} is too large for "
            f"dx={h:g} (lid speed {u_ref:g}). The explicit scheme may diverge "
            f"over long runs; use dt <~ {h / u_ref:.4g} (e.g. {0.5 * h / u_ref:.4g}).",
            stacklevel=2,
        )
    if dt > dt_visc:
        warnings.warn(
            f"Diffusive limit exceeded: dt={dt:g} > Re*dx^2/4 = {dt_visc:.4g}. "
            f"Reduce dt for stability.",
            stacklevel=2,
        )
    return cfl, dt_visc


def simulate_flow(u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve,
                  return_diagnostics=False):

    dx = 1.0/Nxc
    dy = 1.0/Nyc
    x = np.linspace(0 - dx, 1.0 + dx, Nxc+2)
    y = np.linspace(0 - dy, 1.0 + dy, Nyc+2)
    N = tf/dt
    error = []          # pressure-iteration RMSE per timestep
    umax_history = []    # max |u| per timestep (steady-state / divergence monitor)
    diverged_step = None
    t = 0   # initial time

    cfl, dt_visc = _check_stability(dt, dx, dy, Re, Us, Un, Vw, Ve)

    c = np.zeros((Nxc+2, Nyc+2))
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    P = np.zeros((Nxc+2, Nyc+2))
    u_contour = np.zeros((Nxc+1, Nyc+1))
    v_contour = np.zeros((Nxc+1, Nyc+1))
    p_contour = np.zeros((Nxc+1, Nyc+1))
    phi_contour = np.zeros((Nxc+2, Nyc+2))
    phi = np.zeros((Nxc+2, Nyc+2))
    vorticity = np.zeros((Nxc+2, Nyc+2))
    pressure_bc(Nxc, Nyc, c)

    # Pre-compute the red/black interior masks once.
    red, black = checkerboard_masks(Nxc, Nyc)

    # Initialise the plotting grid so a return is always well-defined, even if
    # the very first step diverges.
    [X, Y] = np.meshgrid(x, y)
    STREAMLINE = np.transpose(phi_contour[0:Nxc+2, 0:Nyc+2])

    for n in range(1, int(N+1)):

        print(f'Running iteration number {n} of {int(N)}', end="\r")

        # --- Boundary conditions ---
        velocity_bc(u, v, Us, Un, Vw, Ve, Nxc, Nyc)

        # --- Explicit momentum step ---
        u_tp = navier_stokes_x(u, v, dx, dy, dt, Re)
        v_tp = navier_stokes_y(u, v, dx, dy, dt, Re)

        # --- SOR sweeps for pressure & streamfunction ---
        P_chk = P.copy()
        for _ in range(1, n_iter+1):
            poisson(u_tp, v_tp, dx, dy, dt, P, c, SOR, red, black)
            streamlines(dx, dy, dt, phi, vorticity, SOR, red, black)
            pressure_ghosts(Nxc, Nyc, P)

        rmse = np.sqrt(np.sum((P_chk - P)**2)/(Nxc*Nyc))
        error.append(rmse)

        # --- Vorticity boundary conditions ---
        vorticity_bc(Nxc, Nyc, dx, dy, phi, vorticity)

        # --- Velocity correction (projection) ---
        u[1:Nxc, 1:Nyc+1] = (
            u_tp[1:Nxc, 1:Nyc+1] - (dt/dx)*(P[2:Nxc+1, 1:Nyc+1] - P[1:Nxc, 1:Nyc+1])
        )
        v[1:Nxc+1, 1:Nyc] = (
            v_tp[1:Nxc+1, 1:Nyc] - (dt/dy)*(P[1:Nxc+1, 2:Nyc+1] - P[1:Nxc+1, 1:Nyc])
        )

        # --- Divergence guard ---
        # The explicit scheme can go unstable. Catch it at a sane magnitude
        # (well below float overflow) so we stop *before* NaNs/overflows reach
        # the plots, report when it happened, and return the last finite state.
        # The contour fields below are NOT updated on this step, so the returned
        # state is the previous (still-bounded) one.
        u_ref = max(abs(Un), abs(Us), abs(Vw), abs(Ve), 1.0)
        u_limit = 1.0e3 * u_ref       # lid-cavity velocities are O(1); 1000x = diverged
        umax = np.abs(u).max()
        umax_history.append(umax)
        if not np.isfinite(umax) or umax > u_limit:
            diverged_step = n
            warnings.warn(
                f"Solution diverged at step {n} (t={t + dt:.4g}): max|u|={umax:.3g} "
                f"blew up. dt={dt:g} is too large (convective CFL={cfl:.2f}); "
                f"reduce dt (e.g. {0.5 * min(dx, dy) / u_ref:.4g}). "
                f"Returning the last finite state (step {n - 1}, t={t:.4g}).",
                stacklevel=2,
            )
            break

        # --- Plotting grid ---
        x[0:Nxc+2] = np.linspace(0, Nxc, Nxc+2)
        y[0:Nyc+2] = np.linspace(0, Nyc, Nyc+2)

        # --- Cell-centred contour fields ---
        u_contour[:, :] = 0.5 * (u[0:Nxc+1, 0:Nyc+1] + u[0:Nxc+1, 1:Nyc+2])
        v_contour[:, :] = 0.5 * (v[0:Nxc+1, 0:Nyc+1] + v[1:Nxc+2, 0:Nyc+1])
        p_contour[:, :] = 0.25 * (
            P[0:Nxc+1, 0:Nyc+1] + P[1:Nxc+2, 0:Nyc+1]
            + P[0:Nxc+1, 1:Nyc+2] + P[1:Nxc+2, 1:Nyc+2]
        )
        vorticity[0:Nxc+1, 0:Nyc+1] = (
            (v[1:Nxc+2, 0:Nyc+1] - v[0:Nxc+1, 0:Nyc+1]) / dx
            - (u[0:Nxc+1, 1:Nyc+2] - u[0:Nxc+1, 0:Nyc+1]) / dy
        )
        phi_contour[0:Nxc+1, 0:Nyc+1] = 0.25 * (
            phi[0:Nxc+1, 0:Nyc+1] + phi[1:Nxc+2, 0:Nyc+1]
            + phi[0:Nxc+1, 1:Nyc+2] + phi[1:Nxc+2, 1:Nyc+2]
        )

        [X, Y] = np.meshgrid(x, y)
        STREAMLINE = np.transpose(phi_contour[0:Nxc+2, 0:Nyc+2])

        t += dt

    if return_diagnostics:
        diagnostics = {
            "error": np.asarray(error),
            "umax": np.asarray(umax_history),
            "cfl": cfl,
            "dt_diffusive_limit": dt_visc,
            "diverged_step": diverged_step,
            "diverged_time": None if diverged_step is None else diverged_step * dt,
            "steps_completed": len(umax_history),
        }
        return X, Y, STREAMLINE, u_contour, v_contour, p_contour, vorticity, diagnostics

    return X, Y, STREAMLINE, u_contour, v_contour, p_contour, vorticity
