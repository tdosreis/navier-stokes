"""Lid-driven cavity Navier-Stokes solver.

The coordinate mesh is built once. Contour arrays are refreshed after each
good step so a diverged run can still return the last finite picture. The
vorticity curl stays inside the step because the next streamfunction sweep
reads it.
"""

import warnings

import numpy as np

from boundary_conditions import (
    pressure_bc,
    pressure_ghosts,
    velocity_bc,
    vorticity_bc,
)
from navier_stokes import navier_stokes_x, navier_stokes_y
from solver import (
    checkerboard_masks,
    poisson_sweep,
    source_term,
    stream_sweep,
)


def _check_stability(dt, dx, dy, Re, Us, Un, Vw, Ve):
    """Warn if ``dt`` violates the explicit-scheme stability limits."""

    h = min(dx, dy)
    u_ref = max(abs(Us), abs(Un), abs(Vw), abs(Ve), 1e-30)
    cfl = u_ref * dt / h
    dt_visc = Re * h**2 / 4.0

    if cfl > 1.0:
        warnings.warn(
            f"Convective CFL = {cfl:.2f} > 1: dt={dt:g} is too large for "
            f"dx={h:g} (lid speed {u_ref:g}). The explicit scheme may diverge "
            f"over long runs; use dt <~ {h / u_ref:.4g} (e.g. {0.5 * h / u_ref:.4g}).",
            stacklevel=3,
        )
    if dt > dt_visc:
        warnings.warn(
            f"Diffusive limit exceeded: dt={dt:g} > Re*dx^2/4 = {dt_visc:.4g}. "
            f"Reduce dt for stability.",
            stacklevel=3,
        )
    return cfl, dt_visc


def _curl(Nxc, Nyc, dx, dy, u, v, vorticity):
    """Cell-centred vorticity from the staggered velocity. Feeds the next sweep."""

    vorticity[0:Nxc+1, 0:Nyc+1] = (
        (v[1:Nxc+2, 0:Nyc+1] - v[0:Nxc+1, 0:Nyc+1]) / dx
        - (u[0:Nxc+1, 1:Nyc+2] - u[0:Nxc+1, 0:Nyc+1]) / dy
    )


def _fill_contours(Nxc, Nyc, u, v, P, phi, u_contour, v_contour, p_contour, phi_contour):
    """Cell-centred fields written into the arrays returned to the caller."""

    u_contour[:, :] = 0.5 * (u[0:Nxc+1, 0:Nyc+1] + u[0:Nxc+1, 1:Nyc+2])
    v_contour[:, :] = 0.5 * (v[0:Nxc+1, 0:Nyc+1] + v[1:Nxc+2, 0:Nyc+1])
    p_contour[:, :] = 0.25 * (
        P[0:Nxc+1, 0:Nyc+1] + P[1:Nxc+2, 0:Nyc+1]
        + P[0:Nxc+1, 1:Nyc+2] + P[1:Nxc+2, 1:Nyc+2]
    )
    phi_contour[0:Nxc+1, 0:Nyc+1] = 0.25 * (
        phi[0:Nxc+1, 0:Nyc+1] + phi[1:Nxc+2, 0:Nyc+1]
        + phi[0:Nxc+1, 1:Nyc+2] + phi[1:Nxc+2, 1:Nyc+2]
    )


def _advance(u, v, u_tp, v_tp, P, P_chk, phi, vorticity, c, rhs, neigh,
             red, black, Nxc, Nyc, dx, dy, dt, n_iter, SOR, Re,
             Us, Un, Vw, Ve, u_limit):
    """One timestep. ``u_tp`` / ``v_tp`` / ``rhs`` are reused."""

    velocity_bc(u, v, Us, Un, Vw, Ve, Nxc, Nyc)
    navier_stokes_x(u, v, dx, dy, dt, Re, out=u_tp)
    navier_stokes_y(u, v, dx, dy, dt, Re, out=v_tp)

    P_chk[...] = P
    source_term(u_tp, v_tp, dx, dt, rhs)
    for _ in range(n_iter):
        poisson_sweep(P, c, SOR, red, black, rhs, neigh)
        stream_sweep(phi, vorticity, dx, dy, SOR, red, black, neigh)
        pressure_ghosts(Nxc, Nyc, P)

    rmse = float(np.sqrt(np.sum((P_chk - P)**2) / (Nxc*Nyc)))
    vorticity_bc(Nxc, Nyc, dx, dy, phi, vorticity)

    u[1:Nxc, 1:Nyc+1] = (
        u_tp[1:Nxc, 1:Nyc+1] - (dt/dx)*(P[2:Nxc+1, 1:Nyc+1] - P[1:Nxc, 1:Nyc+1])
    )
    v[1:Nxc+1, 1:Nyc] = (
        v_tp[1:Nxc+1, 1:Nyc] - (dt/dy)*(P[1:Nxc+1, 2:Nyc+1] - P[1:Nxc+1, 1:Nyc])
    )

    umax = float(np.abs(u).max())
    ok = bool(np.isfinite(umax) and umax <= u_limit)
    if ok:
        _curl(Nxc, Nyc, dx, dy, u, v, vorticity)
    return rmse, umax, ok


def simulate_flow(u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve,
                  return_diagnostics=False, callback=None, callback_every=1,
                  verbose=True):
    """Advance the lid-driven cavity flow in time.

    ``callback(n, t, fields)`` is optional and is called every ``callback_every``
    steps (and on the final step). ``fields`` holds the live contour arrays.
    It does not change the solution.

    The passed-in ``u`` and ``v`` are not used as initial data; the cavity
    starts from rest, matching the original script.
    """

    dx = 1.0/Nxc
    dy = 1.0/Nyc
    n_steps = int(tf/dt)
    error = []
    umax_history = []
    diverged_step = None
    t = 0.0

    cfl, dt_visc = _check_stability(dt, dx, dy, Re, Us, Un, Vw, Ve)
    u_ref = max(abs(Un), abs(Us), abs(Vw), abs(Ve), 1.0)
    u_limit = 1.0e3 * u_ref

    c = np.zeros((Nxc+2, Nyc+2))
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    u_tp = np.zeros((Nxc+1, Nyc+2))
    v_tp = np.zeros((Nxc+2, Nyc+1))
    P = np.zeros((Nxc+2, Nyc+2))
    P_chk = np.zeros((Nxc+2, Nyc+2))
    phi = np.zeros((Nxc+2, Nyc+2))
    vorticity = np.zeros((Nxc+2, Nyc+2))
    u_contour = np.zeros((Nxc+1, Nyc+1))
    v_contour = np.zeros((Nxc+1, Nyc+1))
    p_contour = np.zeros((Nxc+1, Nyc+1))
    phi_contour = np.zeros((Nxc+2, Nyc+2))
    pressure_bc(Nxc, Nyc, c)

    red, black = checkerboard_masks(Nxc, Nyc)
    rhs = np.zeros((Nxc+2, Nyc+2))
    neigh = np.zeros((Nxc+2, Nyc+2))

    # Index coordinates, same grid the original loop rebuilt every step.
    # STREAMLINE is a view of phi_contour, so filling contours updates it.
    x = np.linspace(0, Nxc, Nxc+2)
    y = np.linspace(0, Nyc, Nyc+2)
    X, Y = np.meshgrid(x, y)
    STREAMLINE = np.transpose(phi_contour)

    for n in range(1, n_steps+1):
        if verbose:
            print(f"Running iteration number {n} of {n_steps}", end="\r")

        rmse, umax, ok = _advance(
            u, v, u_tp, v_tp, P, P_chk, phi, vorticity, c, rhs, neigh,
            red, black, Nxc, Nyc, dx, dy, dt, n_iter, SOR, Re,
            Us, Un, Vw, Ve, u_limit,
        )

        error.append(rmse)
        umax_history.append(umax)
        if not ok:
            diverged_step = n
            warnings.warn(
                f"Solution diverged at step {n} (t={t + dt:.4g}): max|u|={umax:.3g} "
                f"blew up. dt={dt:g} is too large (convective CFL={cfl:.2f}); "
                f"reduce dt (e.g. {0.5 * min(dx, dy) / u_ref:.4g}). "
                f"Returning the last finite state (step {n - 1}, t={t:.4g}).",
                stacklevel=2,
            )
            break

        t += dt
        _fill_contours(
            Nxc, Nyc, u, v, P, phi, u_contour, v_contour, p_contour, phi_contour
        )
        if callback is not None and (n % callback_every == 0 or n == n_steps):
            callback(n, t, {
                "X": X, "Y": Y, "STREAMLINE": STREAMLINE,
                "u_contour": u_contour, "v_contour": v_contour,
                "p_contour": p_contour, "vorticity": vorticity,
            })

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
