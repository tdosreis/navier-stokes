import numpy as np
import time
import matplotlib.pyplot as plt
from flow_simulation import simulate_flow
from flow_simulation_vectorized import simulate_flow_vectorized
from flow_simulation_optimized import simulate_flow_optimized, simulate_flow_ultra_optimized

def benchmark_implementations(grid_sizes, tf=1.0, n_iter=5):
    """
    Benchmark all implementations across different grid sizes.
    """
    results = {
        'original': [],
        'vectorized': [],
        'optimized': [],
        'ultra_optimized': []
    }
    
    # Test parameters
    Re = 100
    SOR = 1.5
    dt = 0.01
    
    Un = 1.0
    Us = 0.0
    Vw = 0.0
    Ve = 0.0
    
    for Nxc, Nyc in grid_sizes:
        print(f"\nTesting grid size: {Nxc}x{Nyc}")
        print("-" * 40)
        
        u = np.zeros((Nxc+1, Nyc+2))
        v = np.zeros((Nxc+2, Nyc+1))
        
        # Test original implementation
        try:
            print("Testing original implementation...")
            start_time = time.time()
            simulate_flow(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
            orig_time = time.time() - start_time
            results['original'].append(orig_time)
            print(f"Original: {orig_time:.4f}s")
        except Exception as e:
            print(f"Original failed: {e}")
            results['original'].append(np.nan)
        
        # Test vectorized implementation
        try:
            print("Testing vectorized implementation...")
            start_time = time.time()
            simulate_flow_vectorized(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
            vec_time = time.time() - start_time
            results['vectorized'].append(vec_time)
            print(f"Vectorized: {vec_time:.4f}s")
        except Exception as e:
            print(f"Vectorized failed: {e}")
            results['vectorized'].append(np.nan)
        
        # Test optimized implementation
        try:
            print("Testing optimized implementation...")
            start_time = time.time()
            simulate_flow_optimized(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
            opt_time = time.time() - start_time
            results['optimized'].append(opt_time)
            print(f"Optimized: {opt_time:.4f}s")
        except Exception as e:
            print(f"Optimized failed: {e}")
            results['optimized'].append(np.nan)
        
        # Test ultra-optimized implementation
        try:
            print("Testing ultra-optimized implementation...")
            start_time = time.time()
            simulate_flow_ultra_optimized(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
            ultra_time = time.time() - start_time
            results['ultra_optimized'].append(ultra_time)
            print(f"Ultra-optimized: {ultra_time:.4f}s")
        except Exception as e:
            print(f"Ultra-optimized failed: {e}")
            results['ultra_optimized'].append(np.nan)
    
    return results

def plot_performance_results(grid_sizes, results):
    """
    Plot performance comparison results.
    """
    grid_labels = [f"{nx}x{ny}" for nx, ny in grid_sizes]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot execution times
    implementations = ['original', 'vectorized', 'optimized', 'ultra_optimized']
    colors = ['red', 'blue', 'green', 'purple']
    
    for i, impl in enumerate(implementations):
        times = results[impl]
        valid_times = [t for t in times if not np.isnan(t)]
        valid_grids = [j for j, t in enumerate(times) if not np.isnan(t)]
        
        if valid_times:
            ax1.plot(valid_grids, valid_times, 'o-', color=colors[i], label=impl.replace('_', ' ').title(), linewidth=2, markersize=8)
    
    ax1.set_xlabel('Grid Size Index')
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Execution Time Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(range(len(grid_labels)))
    ax1.set_xticklabels(grid_labels, rotation=45)
    
    # Plot speedup relative to original
    for i, impl in enumerate(implementations[1:], 1):  # Skip original
        times = results[impl]
        orig_times = results['original']
        speedups = []
        
        for j, (t, orig_t) in enumerate(zip(times, orig_times)):
            if not np.isnan(t) and not np.isnan(orig_t) and orig_t > 0:
                speedups.append(orig_t / t)
            else:
                speedups.append(np.nan)
        
        valid_speedups = [s for s in speedups if not np.isnan(s)]
        valid_indices = [j for j, s in enumerate(speedups) if not np.isnan(s)]
        
        if valid_speedups:
            ax2.plot(valid_indices, valid_speedups, 'o-', color=colors[i], 
                    label=f"{impl.replace('_', ' ').title()} Speedup", linewidth=2, markersize=8)
    
    ax2.set_xlabel('Grid Size Index')
    ax2.set_ylabel('Speedup Factor')
    ax2.set_title('Speedup Relative to Original')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(range(len(grid_labels)))
    ax2.set_xticklabels(grid_labels, rotation=45)
    
    plt.tight_layout()
    plt.savefig('performance_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def print_summary_table(grid_sizes, results):
    """
    Print a summary table of performance results.
    """
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY")
    print("="*80)
    
    grid_labels = [f"{nx}x{ny}" for nx, ny in grid_sizes]
    implementations = ['original', 'vectorized', 'optimized', 'ultra_optimized']
    
    # Header
    header = f"{'Grid Size':<12}"
    for impl in implementations:
        header += f"{impl.replace('_', ' ').title():<15}"
    print(header)
    print("-" * 80)
    
    # Data rows
    for i, grid_label in enumerate(grid_labels):
        row = f"{grid_label:<12}"
        orig_time = results['original'][i]
        
        for impl in implementations:
            time_val = results[impl][i]
            if np.isnan(time_val):
                row += f"{'N/A':<15}"
            else:
                if impl == 'original':
                    row += f"{time_val:.4f}s{'':<10}"
                else:
                    speedup = orig_time / time_val if not np.isnan(orig_time) and orig_time > 0 else np.nan
                    if np.isnan(speedup):
                        row += f"{time_val:.4f}s{'':<10}"
                    else:
                        row += f"{time_val:.4f}s ({speedup:.1f}x){'':<5}"
        print(row)
    
    # Average speedup
    print("-" * 80)
    print("AVERAGE SPEEDUP:")
    for impl in implementations[1:]:
        speedups = []
        for i in range(len(grid_sizes)):
            orig_time = results['original'][i]
            impl_time = results[impl][i]
            if not np.isnan(orig_time) and not np.isnan(impl_time) and orig_time > 0:
                speedups.append(orig_time / impl_time)
        
        if speedups:
            avg_speedup = np.mean(speedups)
            print(f"{impl.replace('_', ' ').title()}: {avg_speedup:.2f}x")

def main():
    """
    Run comprehensive performance comparison.
    """
    print("="*60)
    print("ADVANCED NAVIER-STOKES PERFORMANCE COMPARISON")
    print("="*60)
    
    # Test grid sizes (smaller for quick testing)
    grid_sizes = [
        (16, 16),
        (32, 32),
        (48, 48),
        (64, 64)
    ]
    
    # Run benchmarks
    results = benchmark_implementations(grid_sizes, tf=0.5, n_iter=3)
    
    # Print summary
    print_summary_table(grid_sizes, results)
    
    # Plot results
    plot_performance_results(grid_sizes, results)
    
    print(f"\n✅ Performance comparison completed!")
    print(f"📊 Results saved to 'performance_comparison.png'")

if __name__ == "__main__":
    main() 