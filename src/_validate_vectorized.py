"""Validation: vectorized backend vs original scalar backend."""
import numpy as np

import navier_stokes as ns_scalar
import navier_stokes_vectorized as ns_vec
from flow_simulation import simulate_flow


def test_momentum_bit_identical():
    rng = np.random.default_rng(0)
    Nxc = Nyc = 24
    dx = dy = 1.0/Nxc
    dt, Re = 0.005, 100.0
    u = rng.standard_normal((Nxc+1, Nyc+2))
    v = rng.standard_normal((Nxc+2, Nyc+1))

    u_tp = np.zeros_like(u)
    for i in range(1, Nxc):
        for j in range(1, Nyc+1):
            ns_scalar.navier_stokes_x(i, j, u, u_tp, v, dx, dy, dt, Re)
    u_tp_vec = ns_vec.navier_stokes_x(u, v, dx, dy, dt, Re)

    v_tp = np.zeros_like(v)
    for i in range(1, Nxc+1):
        for j in range(1, Nyc):
            ns_scalar.navier_stokes_y(i, j, u, v_tp, v, dx, dy, dt, Re)
    v_tp_vec = ns_vec.navier_stokes_y(u, v, dx, dy, dt, Re)

    ex = np.abs(u_tp - u_tp_vec).max()
    ey = np.abs(v_tp - v_tp_vec).max()
    print(f"[momentum] max|du_tp|={ex:.3e}  max|dv_tp|={ey:.3e}")
    assert ex < 1e-12 and ey < 1e-12, "momentum step is NOT bit-identical"


def test_full_pipeline():
    Nxc = Nyc = 32
    Re, tf, SOR, n_iter, dt = 100.0, 0.1, 1.5, 10, 0.01
    Un, Us, Vw, Ve = 1.0, 0.0, 0.0, 0.0
    u0 = np.zeros((Nxc+1, Nyc+2))
    v0 = np.zeros((Nxc+2, Nyc+1))
    args = (Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)

    out_s = simulate_flow(u0, v0, *args, backend="scalar", verbose=False)
    out_v = simulate_flow(u0, v0, *args, backend="vectorized", verbose=False)
    names = ["X", "Y", "STREAMLINE", "u_contour", "v_contour", "p_contour", "vorticity"]
    print("\n[full pipeline] field-by-field comparison (32x32, 10 steps):")
    for nm, a, b in zip(names, out_s, out_v):
        a = np.asarray(a)
        b = np.asarray(b)
        denom = max(np.abs(a).max(), 1e-30)
        print(f"  {nm:11s} max|abs|={np.abs(a-b).max():.3e}  rel={np.abs(a-b).max()/denom:.3e}")


if __name__ == "__main__":
    test_momentum_bit_identical()
    test_full_pipeline()
    print("\nOK")
