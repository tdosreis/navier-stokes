import numpy as np
import time
from flow_simulation_vectorized import simulate_flow_vectorized
from flow_simulation_optimized import simulate_flow_optimized, simulate_flow_ultra_optimized

def test_optimized_consistency():
    """Test that optimized versions produce consistent results."""
    print("Testing optimized versions consistency...")
    
    # Test parameters (small for quick testing)
    Nxc = 16
    Nyc = 16
    Re = 100
    tf = 0.5
    SOR = 1.5
    n_iter = 3
    dt = 0.01
    
    Un = 1.0
    Us = 0.0
    Vw = 0.0
    Ve = 0.0
    
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    
    # Run vectorized version (reference)
    print("Running vectorized version (reference)...")
    start_time = time.time()
    X_ref, Y_ref, STREAMLINE_ref, u_contour_ref, v_contour_ref, p_contour_ref, vorticity_ref = (
        simulate_flow_vectorized(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
    )
    vec_time = time.time() - start_time
    
    # Run optimized version
    print("Running optimized version...")
    start_time = time.time()
    X_opt, Y_opt, STREAMLINE_opt, u_contour_opt, v_contour_opt, p_contour_opt, vorticity_opt = (
        simulate_flow_optimized(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
    )
    opt_time = time.time() - start_time
    
    # Run ultra-optimized version
    print("Running ultra-optimized version...")
    start_time = time.time()
    X_ultra, Y_ultra, STREAMLINE_ultra, u_contour_ultra, v_contour_ultra, p_contour_ultra, vorticity_ultra = (
        simulate_flow_ultra_optimized(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
    )
    ultra_time = time.time() - start_time
    
    # Compare results
    print(f"\nTiming comparison:")
    print(f"Vectorized:      {vec_time:.4f} seconds")
    print(f"Optimized:       {opt_time:.4f} seconds")
    print(f"Ultra-optimized: {ultra_time:.4f} seconds")
    print(f"Optimized speedup:    {vec_time/opt_time:.2f}x")
    print(f"Ultra-optimized speedup: {vec_time/ultra_time:.2f}x")
    
    # Check numerical consistency
    tolerance = 1e-12
    
    # Compare optimized vs vectorized
    opt_consistent = True
    for name, ref, opt in [
        ('u_contour', u_contour_ref, u_contour_opt),
        ('v_contour', v_contour_ref, v_contour_opt),
        ('p_contour', p_contour_ref, p_contour_opt),
        ('streamline', STREAMLINE_ref, STREAMLINE_opt),
        ('vorticity', vorticity_ref, vorticity_opt)
    ]:
        diff = np.max(np.abs(ref - opt))
        if diff > tolerance:
            print(f"❌ {name}: {diff:.2e} (FAIL)")
            opt_consistent = False
        else:
            print(f"✅ {name}: {diff:.2e} (PASS)")
    
    # Compare ultra-optimized vs vectorized
    ultra_consistent = True
    for name, ref, ultra in [
        ('u_contour', u_contour_ref, u_contour_ultra),
        ('v_contour', v_contour_ref, v_contour_ultra),
        ('p_contour', p_contour_ref, p_contour_ultra),
        ('streamline', STREAMLINE_ref, STREAMLINE_ultra),
        ('vorticity', vorticity_ref, vorticity_ultra)
    ]:
        diff = np.max(np.abs(ref - ultra))
        if diff > tolerance:
            print(f"❌ {name}: {diff:.2e} (FAIL)")
            ultra_consistent = False
        else:
            print(f"✅ {name}: {diff:.2e} (PASS)")
    
    print(f"\nConsistency summary:")
    print(f"Optimized:       {'PASS' if opt_consistent else 'FAIL'}")
    print(f"Ultra-optimized: {'PASS' if ultra_consistent else 'FAIL'}")
    
    return opt_consistent and ultra_consistent, vec_time/opt_time, vec_time/ultra_time

if __name__ == "__main__":
    consistent, opt_speedup, ultra_speedup = test_optimized_consistency()
    
    if consistent:
        print(f"\n🎉 All optimized versions are consistent!")
        print(f"🚀 Optimized speedup: {opt_speedup:.1f}x")
        print(f"🚀 Ultra-optimized speedup: {ultra_speedup:.1f}x")
    else:
        print(f"\n❌ Some optimized versions have numerical differences!") 