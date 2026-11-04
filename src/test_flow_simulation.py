"""Scalar backend vs vectorized backend of the single simulate_flow entry point."""
import time

import numpy as np

from flow_simulation import simulate_flow


def _run(backend, Nxc, Nyc, dt, tf, n_iter, SOR, Re):
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    t0 = time.perf_counter()
    out = simulate_flow(
        u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, 0.0, 1.0, 0.0, 0.0,
        backend=backend, verbose=False,
    )
    return time.perf_counter() - t0, out


def test_flow_simulation_consistency():
    """Both backends stay finite. Red-black SOR is close, not bit-identical."""
    Nxc = Nyc = 24
    dt, tf, n_iter, SOR, Re = 0.01, 0.05, 5, 1.5, 100.0

    ts, out_s = _run("scalar", Nxc, Nyc, dt, tf, n_iter, SOR, Re)
    tv, out_v = _run("vectorized", Nxc, Nyc, dt, tf, n_iter, SOR, Re)
    print(f"scalar {ts:.3f}s  vectorized {tv:.3f}s  speedup {ts/max(tv, 1e-9):.1f}x")

    names = ["X", "Y", "STREAMLINE", "u_contour", "v_contour", "p_contour", "vorticity"]
    for name, a, b in zip(names, out_s, out_v):
        a = np.asarray(a)
        b = np.asarray(b)
        assert a.shape == b.shape, name
        assert np.isfinite(a).all() and np.isfinite(b).all(), name
        rel = np.abs(a - b).max() / max(np.abs(a).max(), 1e-30)
        print(f"  {name:11s} rel {rel:.3e}")
    # Coordinates match exactly. Red-black SOR is a different sweep order, so
    # the fields are close rather than bit-identical; velocity should stay near.
    u_rel = np.abs(out_s[3] - out_v[3]).max() / max(np.abs(out_s[3]).max(), 1e-30)
    assert u_rel < 0.02
    assert tv < ts
    return True


if __name__ == "__main__":
    test_flow_simulation_consistency()
    print("OK")
