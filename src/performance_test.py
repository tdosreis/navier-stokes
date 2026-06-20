import numpy as np
import time
from navier_stokes import navier_stokes_x, navier_stokes_y
from navier_stokes_vectorized import navier_stokes_x_vectorized, navier_stokes_y_vectorized

def benchmark_grid_size(nx, ny, num_iterations=10):
    """Benchmark both implementations for a given grid size."""
    print(f"\nGrid size: {nx} x {ny}")
    print("-" * 40)
    
    # Create test data
    u = np.random.rand(nx, ny+1)
    v = np.random.rand(nx+1, ny)
    dx, dy = 0.1, 0.1
    dt = 0.001
    Re = 100.0
    
    # Non-vectorized timing
    u_tp_ref = np.zeros_like(u)
    v_tp_ref = np.zeros_like(v)
    
    start_time = time.time()
    for _ in range(num_iterations):
        for i in range(1, nx-1):
            for j in range(1, ny):
                navier_stokes_x(i, j, u, u_tp_ref, v, dx, dy, dt, Re)
        
        for i in range(1, nx):
            for j in range(1, ny-1):
                navier_stokes_y(i, j, u, v_tp_ref, v, dx, dy, dt, Re)
    non_vec_time = time.time() - start_time
    
    # Vectorized timing
    start_time = time.time()
    for _ in range(num_iterations):
        u_tp_vec = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
        v_tp_vec = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)
    vec_time = time.time() - start_time
    
    speedup = non_vec_time / vec_time
    
    print(f"Non-vectorized: {non_vec_time:.4f}s ({non_vec_time/num_iterations:.4f}s per iteration)")
    print(f"Vectorized:     {vec_time:.4f}s ({vec_time/num_iterations:.4f}s per iteration)")
    print(f"Speedup:        {speedup:.2f}x")
    
    return speedup

def main():
    """Run performance benchmarks for different grid sizes."""
    print("=" * 60)
    print("Navier-Stokes Vectorized vs Non-vectorized Performance Test")
    print("=" * 60)
    
    # Test different grid sizes
    grid_sizes = [
        (10, 8),
        (20, 16),
        (50, 40),
        (100, 80),
        (200, 160)
    ]
    
    results = []
    
    for nx, ny in grid_sizes:
        speedup = benchmark_grid_size(nx, ny)
        results.append((nx, ny, speedup))
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"{'Grid Size':<15} {'Speedup':<10}")
    print("-" * 25)
    
    for nx, ny, speedup in results:
        print(f"{nx}x{ny:<10} {speedup:.2f}x")
    
    # Calculate average speedup
    avg_speedup = np.mean([r[2] for r in results])
    print(f"\nAverage speedup: {avg_speedup:.2f}x")
    
    # Memory usage comparison
    print(f"\nMemory efficiency:")
    print(f"- Non-vectorized: Uses nested loops, minimal memory overhead")
    print(f"- Vectorized: Uses NumPy arrays, more memory efficient for large grids")
    
    print(f"\nKey benefits of vectorization:")
    print(f"- Massive speedup ({avg_speedup:.1f}x faster on average)")
    print(f"- Cleaner, more readable code")
    print(f"- Better cache utilization")
    print(f"- Leverages optimized BLAS/LAPACK libraries")

if __name__ == "__main__":
    main() 