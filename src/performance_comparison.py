import numpy as np
import time
import matplotlib.pyplot as plt

# Import original functions
from navier_stokes import navier_stokes_x, navier_stokes_y
from solver import poisson, streamlines
from boundary_conditions import (
    velocity_bc, pressure_bc, pressure_ghosts, vorticity_bc
)
from flow_simulation import simulate_flow

# Import vectorized functions
from navier_stokes_vectorized import navier_stokes_x_vectorized, navier_stokes_y_vectorized
from solver_vectorized import poisson_vectorized, streamlines_vectorized
from boundary_conditions_vectorized import (
    velocity_bc_vectorized, pressure_bc_vectorized, 
    pressure_ghosts_vectorized, vorticity_bc_vectorized
)
from flow_simulation_vectorized import simulate_flow_vectorized


def benchmark_navier_stokes():
    """Benchmark the Navier-Stokes functions."""
    print("Benchmarking Navier-Stokes functions...")
    
    # Test parameters
    Nxc, Nyc = 50, 50
    dx, dy = 1.0/Nxc, 1.0/Nyc
    dt = 0.001
    Re = 100.0
    
    # Initialize arrays
    u = np.random.rand(Nxc+1, Nyc+2)
    v = np.random.rand(Nxc+2, Nyc+1)
    u_tp = np.zeros((Nxc+1, Nyc+2))
    v_tp = np.zeros((Nxc+2, Nyc+1))
    
    # Benchmark original version
    start_time = time.time()
    for i in range(1, Nxc):
        for j in range(1, Nyc+1):
            u_tp = navier_stokes_x(i, j, u, u_tp, v, dx, dy, dt, Re)
    
    for i in range(1, Nxc+1):
        for j in range(1, Nyc):
            v_tp = navier_stokes_y(i, j, u, v_tp, v, dx, dy, dt, Re)
    
    original_time = time.time() - start_time
    
    # Benchmark vectorized version
    start_time = time.time()
    u_tp_vec = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
    v_tp_vec = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)
    vectorized_time = time.time() - start_time
    
    print(f"Original Navier-Stokes time: {original_time:.4f} seconds")
    print(f"Vectorized Navier-Stokes time: {vectorized_time:.4f} seconds")
    print(f"Speedup: {original_time/vectorized_time:.2f}x")
    
    # Check if results are similar
    u_diff = np.max(np.abs(u_tp - u_tp_vec))
    v_diff = np.max(np.abs(v_tp - v_tp_vec))
    print(f"Maximum difference in u: {u_diff:.2e}")
    print(f"Maximum difference in v: {v_diff:.2e}")
    print()


def benchmark_solver():
    """Benchmark the solver functions."""
    print("Benchmarking solver functions...")
    
    # Test parameters
    Nxc, Nyc = 50, 50
    dx, dy = 1.0/Nxc, 1.0/Nyc
    dt = 0.001
    SOR = 1.5
    
    # Initialize arrays
    u_tp = np.random.rand(Nxc+1, Nyc+2)
    v_tp = np.random.rand(Nxc+2, Nyc+1)
    P = np.random.rand(Nxc+2, Nyc+2)
    c = np.random.rand(Nxc+2, Nyc+2)
    phi = np.random.rand(Nxc+2, Nyc+2)
    vorticity = np.random.rand(Nxc+2, Nyc+2)
    
    # Benchmark original version
    start_time = time.time()
    for i in range(1, Nxc+1):
        for j in range(1, Nyc+1):
            P = poisson(i, j, u_tp, v_tp, dx, dy, dt, P, c, SOR)
            phi = streamlines(i, j, dx, dy, dt, phi, vorticity, SOR)
    
    original_time = time.time() - start_time
    
    # Benchmark vectorized version
    start_time = time.time()
    P_vec = poisson_vectorized(u_tp, v_tp, dx, dy, dt, P.copy(), c, SOR)
    phi_vec = streamlines_vectorized(dx, dy, dt, phi.copy(), vorticity, SOR)
    vectorized_time = time.time() - start_time
    
    print(f"Original solver time: {original_time:.4f} seconds")
    print(f"Vectorized solver time: {vectorized_time:.4f} seconds")
    print(f"Speedup: {original_time/vectorized_time:.2f}x")
    print()


def benchmark_boundary_conditions():
    """Benchmark the boundary condition functions."""
    print("Benchmarking boundary condition functions...")
    
    # Test parameters
    Nxc, Nyc = 50, 50
    dx, dy = 1.0/Nxc, 1.0/Nyc
    Us, Un, Vw, Ve = 1.0, 0.0, 0.0, 0.0
    
    # Initialize arrays
    u = np.random.rand(Nxc+1, Nyc+2)
    v = np.random.rand(Nxc+2, Nyc+1)
    c = np.zeros((Nxc+2, Nyc+2))
    P = np.random.rand(Nxc+2, Nyc+2)
    phi = np.random.rand(Nxc+2, Nyc+2)
    vorticity = np.random.rand(Nxc+2, Nyc+2)
    
    # Benchmark original version
    start_time = time.time()
    pressure_bc(Nxc, Nyc, c)
    for i in range(0, Nxc+1):
        for j in range(0, Nyc+1):
            velocity_bc(i, j, u, v, Us, Un, Vw, Ve, Nxc, Nyc)
    pressure_ghosts(Nxc, Nyc, P)
    for i in range(1, Nxc+1):
        for j in range(1, Nyc+1):
            vorticity_bc(i, j, Nxc, Nyc, dx, dy, phi, vorticity)
    
    original_time = time.time() - start_time
    
    # Benchmark vectorized version
    start_time = time.time()
    pressure_bc_vectorized(Nxc, Nyc, c.copy())
    velocity_bc_vectorized(u.copy(), v.copy(), Us, Un, Vw, Ve, Nxc, Nyc)
    pressure_ghosts_vectorized(Nxc, Nyc, P.copy())
    vorticity_bc_vectorized(Nxc, Nyc, dx, dy, phi.copy(), vorticity.copy())
    vectorized_time = time.time() - start_time
    
    print(f"Original boundary conditions time: {original_time:.4f} seconds")
    print(f"Vectorized boundary conditions time: {vectorized_time:.4f} seconds")
    print(f"Speedup: {original_time/vectorized_time:.2f}x")
    print()


def benchmark_full_simulation():
    """Benchmark the full simulation."""
    print("Benchmarking full simulation...")
    
    # Test parameters (smaller for quick test)
    Nxc, Nyc = 20, 20
    dt = 0.01
    tf = 0.1
    n_iter = 5
    SOR = 1.5
    Re = 100.0
    Us, Un, Vw, Ve = 1.0, 0.0, 0.0, 0.0
    
    # Initialize arrays
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    
    # Benchmark original version
    start_time = time.time()
    try:
        X_orig, Y_orig, STREAMLINE_orig, u_contour_orig, v_contour_orig, p_contour_orig, vorticity_orig = simulate_flow(
            u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve
        )
        original_time = time.time() - start_time
        print(f"Original simulation time: {original_time:.4f} seconds")
    except Exception as e:
        print(f"Original simulation failed: {e}")
        original_time = float('inf')
    
    # Benchmark vectorized version
    start_time = time.time()
    try:
        X_vec, Y_vec, STREAMLINE_vec, u_contour_vec, v_contour_vec, p_contour_vec, vorticity_vec = simulate_flow_vectorized(
            u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve
        )
        vectorized_time = time.time() - start_time
        print(f"Vectorized simulation time: {vectorized_time:.4f} seconds")
        
        if original_time != float('inf'):
            print(f"Speedup: {original_time/vectorized_time:.2f}x")
            
            # Check if results are similar
            u_diff = np.max(np.abs(u_contour_orig - u_contour_vec))
            v_diff = np.max(np.abs(v_contour_orig - v_contour_vec))
            p_diff = np.max(np.abs(p_contour_orig - p_contour_vec))
            print(f"Maximum difference in u_contour: {u_diff:.2e}")
            print(f"Maximum difference in v_contour: {v_diff:.2e}")
            print(f"Maximum difference in p_contour: {p_diff:.2e}")
            
    except Exception as e:
        print(f"Vectorized simulation failed: {e}")
    
    print()


if __name__ == "__main__":
    print("Performance Comparison: Original vs Vectorized Navier-Stokes Solver")
    print("=" * 70)
    
    benchmark_navier_stokes()
    benchmark_solver()
    benchmark_boundary_conditions()
    benchmark_full_simulation()
    
    print("Benchmarking complete!") 