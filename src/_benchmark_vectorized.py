"""Benchmark: scalar backend vs vectorized backend of simulate_flow."""
import time
import numpy as np

from flow_simulation import simulate_flow


def run(backend, Nxc, Nyc, dt, tf, n_iter, SOR, Re):
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    t0 = time.perf_counter()
    simulate_flow(
        u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, 0.0, 1.0, 0.0, 0.0,
        backend=backend, verbose=False,
    )
    return time.perf_counter() - t0


if __name__ == "__main__":
    Re, SOR, n_iter, dt, tf = 1000.0, 1.5, 10, 0.01, 0.1   # 10 timesteps
    print(f"{'grid':>10} {'scalar (s)':>12} {'vector (s)':>12} {'speedup':>9}")
    for N in (16, 32, 48, 64):
        ts = run("scalar", N, N, dt, tf, n_iter, SOR, Re)
        tv = run("vectorized", N, N, dt, tf, n_iter, SOR, Re)
        print(f"{N}x{N:<6} {ts:12.3f} {tv:12.3f} {ts/tv:8.1f}x")
