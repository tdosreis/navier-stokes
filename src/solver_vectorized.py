"""Vectorized pressure / streamfunction SOR solvers.

The scalar ``poisson`` / ``streamlines`` routines sweep the grid cell by cell
and update in place (lexicographic Gauss-Seidel with over-relaxation), which is
inherently sequential. The vectorized equivalent uses **red-black (checkerboard)
ordering**: all "red" cells (i+j even) are updated together using their black
neighbours, then all "black" cells are updated using the freshly updated red
ones. This keeps the Gauss-Seidel character (each colour sees updated
neighbours) while letting NumPy process a whole colour at once, and it converges
to the same solution as the lexicographic sweep.
"""

import numpy as np


def checkerboard_masks(Nxc, Nyc):
    """Return boolean (red, black) interior masks for a (Nxc+2, Nyc+2) field."""

    ii, jj = np.meshgrid(np.arange(Nxc+2), np.arange(Nyc+2), indexing="ij")
    interior = (ii >= 1) & (ii <= Nxc) & (jj >= 1) & (jj <= Nyc)
    red = interior & (((ii + jj) % 2) == 0)
    black = interior & (((ii + jj) % 2) == 1)
    return red, black


def _laplacian_sum(P):
    """Sum of the four orthogonal neighbours, full-size (zeros on the border)."""

    out = np.zeros_like(P)
    out[1:-1, 1:-1] = (
        P[2:, 1:-1] + P[:-2, 1:-1] + P[1:-1, 2:] + P[1:-1, :-2]
    )
    return out


def poisson(u_tp, v_tp, dx, dy, dt, P, c, SOR, red, black):
    """One red-black SOR sweep of the pressure-Poisson equation (updates ``P``)."""

    # Mass-conservation source term (constant across the inner sweeps).
    # Interior pressure cells are i = 1..Nxc, j = 1..Nyc. The staggered u_tp
    # (Nxc+1, Nyc+2) and v_tp (Nxc+2, Nyc+1) arrays need their own slices.
    Nxc, Nyc = P.shape[0] - 2, P.shape[1] - 2
    rhs = np.zeros_like(P)
    rhs[1:Nxc+1, 1:Nyc+1] = -(dx/dt) * (
        u_tp[1:Nxc+1, 1:Nyc+1] - u_tp[0:Nxc, 1:Nyc+1]
        + v_tp[1:Nxc+1, 1:Nyc+1] - v_tp[1:Nxc+1, 0:Nyc]
    )

    for colour in (red, black):
        neigh = _laplacian_sum(P)
        P[colour] = (
            SOR * c[colour] * (neigh[colour] + rhs[colour])
            + (1.0 - SOR) * P[colour]
        )

    return P


def streamlines(dx, dy, dt, phi, vorticity, SOR, red, black):
    """One red-black SOR sweep of the streamfunction equation (updates ``phi``)."""

    for colour in (red, black):
        neigh = _laplacian_sum(phi)
        phi[colour] = (
            0.25 * SOR * (neigh[colour] + (dx*dy) * vorticity[colour])
            + (1.0 - SOR) * phi[colour]
        )

    return phi
