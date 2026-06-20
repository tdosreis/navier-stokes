"""Vectorized boundary conditions for the lid-driven cavity solver.

These mirror the scalar routines in ``boundary_conditions.py`` but apply each
condition to whole array slices at once instead of looping over every cell.
"""


def pressure_bc(Nxc, Nyc, c):
    """Build the SOR coefficient field ``c`` (1/n_neighbours per cell)."""

    c[2:Nxc, 2:Nyc] = 1.0/4.0
    c[1, 1] = 1.0/2.0
    c[Nxc, 1] = 1.0/2.0
    c[1, Nyc] = 1.0/2.0
    c[Nxc, Nyc] = 1.0/2.0
    c[2:Nxc, 1] = 1.0/3.0
    c[1, 2:Nyc] = 1.0/3.0
    c[Nxc, 2:Nyc] = 1.0/3.0
    c[2:Nxc, Nyc] = 1.0/3.0

    return c


def velocity_bc(u, v, Us, Un, Vw, Ve, Nxc, Nyc):
    """Apply ghost-cell velocity boundary conditions to ``u`` and ``v``.

    This reproduces the exact net result of the scalar nested-loop routine,
    including its corner behaviour: the left/right wall columns of ``u`` are
    fully zeroed (so the two top-corner ghosts are 0, not ``2*Un``), and the
    right ghost column of ``v`` is zeroed at the bottom/top corners.
    """

    # No-penetration walls: zero the left/right ghost columns of u and the
    # bottom/top ghost rows of v first, so the ghost expressions below pick up
    # the zeroed neighbours exactly as the scalar loop does.
    u[0, 1:Nyc+2] = 0.0
    u[Nxc, 1:Nyc+2] = 0.0
    v[1:Nxc+1, 0] = 0.0
    v[1:Nxc+1, Nyc] = 0.0

    # Bottom / top ghosts for u over the interior columns (moving lid = Un).
    u[1:Nxc, 0] = 2.0 * Us - u[1:Nxc, 1]
    u[1:Nxc, Nyc+1] = 2.0 * Un - u[1:Nxc, Nyc]
    u[0, 0] = 2.0 * Us - u[0, 1]
    u[Nxc, 0] = 2.0 * Us - u[Nxc, 1]

    # Left / right ghost columns for v.
    v[0, :] = 2.0 * Vw - v[1, :]
    v[Nxc+1, :] = 2.0 * Ve - v[Nxc, :]
    v[Nxc+1, 0] = 0.0
    v[Nxc+1, Nyc] = 0.0

    return v, u


def pressure_ghosts(Nxc, Nyc, P):
    """Zero the pressure ghost cells around the domain."""

    P[0:Nxc+2, 0] = 0
    P[0:Nxc+2, Nyc+1] = 0
    P[0, 0:Nyc+2] = 0
    P[Nxc+1, 0:Nyc+2] = 0

    return P


def vorticity_bc(Nxc, Nyc, dx, dy, phi, vorticity):
    """Set the vorticity on the four domain boundaries from the streamfunction."""

    vorticity[1:Nxc+1, 0] = -2.0 * phi[1:Nxc+1, 1] / (dx*dy)
    vorticity[1:Nxc+1, Nyc+1] = -2.0 * phi[1:Nxc+1, Nyc] / (dx*dy) - 2.0/dx
    vorticity[0, 1:Nyc+1] = -2.0 * phi[1, 1:Nyc+1] / (dx*dy)
    vorticity[Nxc+1, 1:Nyc+1] = -2.0 * phi[Nxc, 1:Nyc+1] / (dx*dy)

    return vorticity
