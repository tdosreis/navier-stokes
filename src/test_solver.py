"""Checks for the single-entry-point solver."""
import warnings

import numpy as np

from flow_simulation import simulate_flow
from navier_stokes import navier_stokes_x, navier_stokes_y
from navier_stokes_vectorized import navier_stokes_x as nsx
from navier_stokes_vectorized import navier_stokes_y as nsy


def test_momentum_buffer_matches_scalar():
    rng = np.random.default_rng(1)
    Nxc = Nyc = 16
    u = rng.standard_normal((Nxc+1, Nyc+2))
    v = rng.standard_normal((Nxc+2, Nyc+1))
    dx = dy = 1.0/Nxc
    dt, Re = 0.002, 100.0

    u_ref = np.zeros_like(u)
    for i in range(1, Nxc):
        for j in range(1, Nyc+1):
            navier_stokes_x(i, j, u, u_ref, v, dx, dy, dt, Re)
    buf = np.zeros_like(u)
    u_out = nsx(u, v, dx, dy, dt, Re, out=buf)
    assert u_out is buf
    assert np.abs(u_ref - u_out).max() < 1e-12

    v_ref = np.zeros_like(v)
    for i in range(1, Nxc+1):
        for j in range(1, Nyc):
            navier_stokes_y(i, j, u, v_ref, v, dx, dy, dt, Re)
    buf = np.zeros_like(v)
    v_out = nsy(u, v, dx, dy, dt, Re, out=buf)
    assert np.abs(v_ref - v_out).max() < 1e-12


def test_callback_does_not_change_result():
    N, dt, tf = 12, 0.01, 0.04
    args = (N, N, dt, tf, 4, 1.5, 100.0, 0.0, 1.0, 0.0, 0.0)
    seen = []

    def cb(n, t, fields):
        seen.append((n, t, float(np.max(np.abs(fields["u_contour"])))))

    u = np.zeros((N+1, N+2))
    v = np.zeros((N+2, N+1))
    plain = simulate_flow(u, v, *args, verbose=False)
    watched = simulate_flow(u, v, *args, callback=cb, callback_every=2, verbose=False)
    for a, b in zip(plain, watched):
        assert np.allclose(a, b)
    assert seen and seen[-1][0] == int(tf/dt)


def test_diagnostics_length():
    N, dt, tf = 8, 0.01, 0.03
    u = np.zeros((N+1, N+2))
    v = np.zeros((N+2, N+1))
    *_, diag = simulate_flow(
        u, v, N, N, dt, tf, 3, 1.5, 100.0, 0.0, 1.0, 0.0, 0.0,
        return_diagnostics=True, verbose=False,
    )
    assert diag["diverged_step"] is None
    assert len(diag["error"]) == int(tf/dt)
    assert len(diag["umax"]) == int(tf/dt)
    assert diag["cfl"] < 1.0


def test_both_backends_stop_when_velocity_blows_up():
    """A huge but finite velocity is rejected, and the picture stays at the last good step."""
    N, dt, tf, Re = 16, 0.1, 0.5, 10.0
    u_limit = 1.0e3  # lid speed is 1, so the shared cap is 1e3 * u_ref
    for backend in ("scalar", "vectorized"):
        u = np.zeros((N+1, N+2))
        v = np.zeros((N+2, N+1))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            *fields, diag = simulate_flow(
                u, v, N, N, dt, tf, 4, 1.5, Re, 0.0, 1.0, 0.0, 0.0,
                return_diagnostics=True, verbose=False, backend=backend,
            )
        assert diag["diverged_step"] is not None, backend
        assert diag["umax"][-1] > u_limit, backend
        u_contour = fields[3]
        assert np.isfinite(u_contour).all(), backend
        assert np.abs(u_contour).max() <= u_limit, backend


if __name__ == "__main__":
    test_momentum_buffer_matches_scalar()
    test_callback_does_not_change_result()
    test_diagnostics_length()
    test_both_backends_stop_when_velocity_blows_up()
    print("OK")
