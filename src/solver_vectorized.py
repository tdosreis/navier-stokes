"""Vectorized pressure / streamfunction SOR solvers.

The scalar ``poisson`` / ``streamlines`` routines sweep the grid cell by cell
and update in place (lexicographic Gauss-Seidel with over-relaxation), which is
inherently sequential. The vectorized equivalent uses **red-black (checkerboard)
ordering**: all "red" cells (i+j even) are updated together using their black
neighbours, then all "black" cells are updated using the freshly updated red
ones. This keeps the Gauss-Seidel character (each colour sees updated
neighbours) while letting NumPy process a whole colour at once, and it converges
to the same solution as the lexicographic sweep.

Work arrays (``rhs``, ``neigh``) are filled in place so a timestep does not
allocate a new grid on every sweep.
"""

import numpy as np


def checkerboard_masks(Nxc, Nyc):
    """Return boolean (red, black) interior masks for a (Nxc+2, Nyc+2) field."""

    ii, jj = np.meshgrid(np.arange(Nxc+2), np.arange(Nyc+2), indexing="ij")
    interior = (ii >= 1) & (ii <= Nxc) & (jj >= 1) & (jj <= Nyc)
    red = interior & (((ii + jj) % 2) == 0)
    black = interior & (((ii + jj) % 2) == 1)
    return red, black


def _laplacian_sum(field, out):
    """Sum of the four orthogonal neighbours. ``out`` is written in place."""

    out[1:-1, 1:-1] = (
        field[2:, 1:-1] + field[:-2, 1:-1] + field[1:-1, 2:] + field[1:-1, :-2]
    )
    return out


def source_term(u_tp, v_tp, dx, dt, rhs):
    """Mass-conservation source for the pressure Poisson equation.

    Constant for every inner SOR sweep of one timestep. Interior pressure cells
    are i = 1..Nxc, j = 1..Nyc.
    """

    Nxc, Nyc = rhs.shape[0] - 2, rhs.shape[1] - 2
    rhs[1:Nxc+1, 1:Nyc+1] = -(dx/dt) * (
        u_tp[1:Nxc+1, 1:Nyc+1] - u_tp[0:Nxc, 1:Nyc+1]
        + v_tp[1:Nxc+1, 1:Nyc+1] - v_tp[1:Nxc+1, 0:Nyc]
    )
    return rhs


def poisson_sweep(P, c, SOR, red, black, rhs, neigh):
    """One red-black SOR sweep of the pressure equation. Updates ``P``."""

    for colour in (red, black):
        _laplacian_sum(P, neigh)
        P[colour] = (
            SOR * c[colour] * (neigh[colour] + rhs[colour])
            + (1.0 - SOR) * P[colour]
        )
    return P


def stream_sweep(phi, vorticity, dx, dy, SOR, red, black, neigh):
    """One red-black SOR sweep of the streamfunction equation. Updates ``phi``."""

    for colour in (red, black):
        _laplacian_sum(phi, neigh)
        phi[colour] = (
            0.25 * SOR * (neigh[colour] + (dx*dy) * vorticity[colour])
            + (1.0 - SOR) * phi[colour]
        )
    return phi
