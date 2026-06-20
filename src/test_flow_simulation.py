import numpy as np
import time
from flow_simulation import simulate_flow
from flow_simulation_vectorized import simulate_flow_vectorized

def test_flow_simulation_consistency():
    """Test that vectorized and non-vectorized flow simulations produce identical results."""
    print("Testing flow simulation consistency...")
    
    # Test parameters (smaller for quick testing)
    Nxc = 32
    Nyc = 32
    Re = 100
    tf = 1.0  # Shorter time for testing
    SOR = 1.5
    n_iter = 5
    dt = 0.01
    
    Un = 1.0
    Us = 0.0
    Vw = 0.0
    Ve = 0.0
    
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    
    # Run non-vectorized version
    print("Running non-vectorized simulation...")
    start_time = time.time()
    X_ref, Y_ref, STREAMLINE_ref, u_contour_ref, v_contour_ref, p_contour_ref, vorticity_ref = (
        simulate_flow(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
    )
    non_vec_time = time.time() - start_time
    
    # Run vectorized version
    print("\nRunning vectorized simulation...")
    start_time = time.time()
    X_vec, Y_vec, STREAMLINE_vec, u_contour_vec, v_contour_vec, p_contour_vec, vorticity_vec = (
        simulate_flow_vectorized(u.copy(), v.copy(), Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
    )
    vec_time = time.time() - start_time
    
    # Compare results
    print(f"\nTiming comparison:")
    print(f"Non-vectorized: {non_vec_time:.4f} seconds")
    print(f"Vectorized:     {vec_time:.4f} seconds")
    print(f"Speedup:        {non_vec_time/vec_time:.2f}x")
    
    # Check numerical consistency with more realistic tolerance
    # Use relative tolerance based on the magnitude of the values
    tolerance = 1e-12
    
    u_diff = np.max(np.abs(u_contour_ref - u_contour_vec))
    v_diff = np.max(np.abs(v_contour_ref - v_contour_vec))
    p_diff = np.max(np.abs(p_contour_ref - p_contour_vec))
    stream_diff = np.max(np.abs(STREAMLINE_ref - STREAMLINE_vec))
    vort_diff = np.max(np.abs(vorticity_ref - vorticity_vec))
    
    # Calculate relative differences
    u_rel_diff = u_diff / (np.max(np.abs(u_contour_ref)) + 1e-16)
    v_rel_diff = v_diff / (np.max(np.abs(v_contour_ref)) + 1e-16)
    p_rel_diff = p_diff / (np.max(np.abs(p_contour_ref)) + 1e-16)
    stream_rel_diff = stream_diff / (np.max(np.abs(STREAMLINE_ref)) + 1e-16)
    vort_rel_diff = vort_diff / (np.max(np.abs(vorticity_ref)) + 1e-16)
    
    print(f"\nAbsolute differences:")
    print(f"u_contour:      {u_diff:.2e}")
    print(f"v_contour:      {v_diff:.2e}")
    print(f"p_contour:      {p_diff:.2e}")
    print(f"streamline:     {stream_diff:.2e}")
    print(f"vorticity:      {vort_diff:.2e}")
    
    print(f"\nRelative differences:")
    print(f"u_contour:      {u_rel_diff:.2e}")
    print(f"v_contour:      {v_rel_diff:.2e}")
    print(f"p_contour:      {p_rel_diff:.2e}")
    print(f"streamline:     {stream_rel_diff:.2e}")
    print(f"vorticity:      {vort_rel_diff:.2e}")
    
    # Check if all differences are within tolerance
    all_consistent = (
        u_diff < tolerance and 
        v_diff < tolerance and 
        p_diff < tolerance and 
        stream_diff < tolerance and 
        vort_diff < tolerance
    )
    
    # Also check relative differences
    rel_tolerance = 1e-10
    all_rel_consistent = (
        u_rel_diff < rel_tolerance and 
        v_rel_diff < rel_tolerance and 
        p_rel_diff < rel_tolerance and 
        stream_rel_diff < rel_tolerance and 
        vort_rel_diff < rel_tolerance
    )
    
    print(f"\nAbsolute consistency check: {'PASS' if all_consistent else 'FAIL'}")
    print(f"Relative consistency check: {'PASS' if all_rel_consistent else 'FAIL'}")
    
    return all_consistent or all_rel_consistent, non_vec_time/vec_time

def test_larger_simulation():
    """Test performance on a larger grid."""
    print("\n" + "="*50)
    print("Testing larger simulation performance...")
    print("="*50)
    
    # Larger test parameters
    Nxc = 64
    Nyc = 64
    Re = 1000
    tf = 2.0
    SOR = 1.5
    n_iter = 10
    dt = 0.01
    
    Un = 1.0
    Us = 0.0
    Vw = 0.0
    Ve = 0.0
    
    u = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    
    # Run vectorized version only (for performance testing)
    print(f"Running vectorized simulation on {Nxc}x{Nyc} grid...")
    start_time = time.time()
    X, Y, STREAMLINE, u_contour, v_contour, p_contour, vorticity = (
        simulate_flow_vectorized(u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve)
    )
    vec_time = time.time() - start_time
    
    print(f"Vectorized simulation completed in {vec_time:.2f} seconds")
    print(f"Grid size: {Nxc}x{Nyc}")
    print(f"Time steps: {int(tf/dt)}")
    print(f"Total iterations: {int(tf/dt) * n_iter}")
    
    return vec_time

if __name__ == "__main__":
    # Test consistency
    consistent, speedup = test_flow_simulation_consistency()
    
    if consistent:
        print(f"\n✅ Vectorized flow simulation is consistent with original!")
        print(f"🚀 Achieved {speedup:.1f}x speedup")
        
        # Test larger simulation
        test_larger_simulation()
    else:
        print(f"\n❌ Vectorized flow simulation has numerical differences!")
        print("Note: Small differences may be due to floating-point arithmetic order differences.") 