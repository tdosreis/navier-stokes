import numpy as np
import time
from navier_stokes import navier_stokes_x, navier_stokes_y
from navier_stokes_vectorized import navier_stokes_x_vectorized, navier_stokes_y_vectorized

def create_test_data(nx=10, ny=8):
    """Create test velocity fields with appropriate shapes."""
    # u shape: (Nxc+1, Nyc+2) = (nx, ny+1)
    # v shape: (Nxc+2, Nyc+1) = (nx+1, ny)
    u = np.random.rand(nx, ny+1)
    v = np.random.rand(nx+1, ny)
    
    # Create output arrays
    u_tp = np.zeros_like(u)
    v_tp = np.zeros_like(v)
    
    return u, v, u_tp, v_tp

def test_consistency():
    """Test that vectorized and non-vectorized versions produce identical results."""
    print("Testing consistency between vectorized and non-vectorized versions...")
    
    # Test parameters
    nx, ny = 10, 8
    dx, dy = 0.1, 0.1
    dt = 0.001
    Re = 100.0
    
    # Create test data
    u, v, u_tp_ref, v_tp_ref = create_test_data(nx, ny)
    u_tp_vec = np.zeros_like(u)
    v_tp_vec = np.zeros_like(v)
    
    # Run non-vectorized version
    for i in range(1, nx-1):
        for j in range(1, ny):
            navier_stokes_x(i, j, u, u_tp_ref, v, dx, dy, dt, Re)
    
    for i in range(1, nx):
        for j in range(1, ny-1):
            navier_stokes_y(i, j, u, v_tp_ref, v, dx, dy, dt, Re)
    
    # Run vectorized version
    u_tp_vec = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
    v_tp_vec = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)
    
    # Compare results
    u_diff = np.abs(u_tp_ref - u_tp_vec)
    v_diff = np.abs(v_tp_ref - v_tp_vec)
    
    u_max_diff = np.max(u_diff)
    v_max_diff = np.max(v_diff)
    
    print(f"Maximum difference in u: {u_max_diff:.2e}")
    print(f"Maximum difference in v: {v_max_diff:.2e}")
    
    # Check if differences are within numerical precision
    tolerance = 1e-14
    u_consistent = u_max_diff < tolerance
    v_consistent = v_max_diff < tolerance
    
    print(f"u consistency: {'PASS' if u_consistent else 'FAIL'}")
    print(f"v consistency: {'PASS' if v_consistent else 'FAIL'}")
    
    return u_consistent and v_consistent

def test_performance():
    """Test performance improvement of vectorized version."""
    print("\nTesting performance...")
    
    # Test parameters
    nx, ny = 100, 80
    dx, dy = 0.1, 0.1
    dt = 0.001
    Re = 100.0
    
    # Create test data
    u, v, _, _ = create_test_data(nx, ny)
    
    # Time non-vectorized version
    start_time = time.time()
    u_tp_ref = np.zeros_like(u)
    v_tp_ref = np.zeros_like(v)
    
    for i in range(1, nx-1):
        for j in range(1, ny):
            navier_stokes_x(i, j, u, u_tp_ref, v, dx, dy, dt, Re)
    
    for i in range(1, nx):
        for j in range(1, ny-1):
            navier_stokes_y(i, j, u, v_tp_ref, v, dx, dy, dt, Re)
    
    non_vec_time = time.time() - start_time
    
    # Time vectorized version
    start_time = time.time()
    u_tp_vec = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
    v_tp_vec = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)
    vec_time = time.time() - start_time
    
    speedup = non_vec_time / vec_time
    
    print(f"Non-vectorized time: {non_vec_time:.4f} seconds")
    print(f"Vectorized time: {vec_time:.4f} seconds")
    print(f"Speedup: {speedup:.2f}x")
    
    return speedup > 1.0

def test_edge_cases():
    """Test edge cases and boundary conditions."""
    print("\nTesting edge cases...")
    
    # Test with small grid
    nx, ny = 3, 3
    dx, dy = 0.1, 0.1
    dt = 0.001
    Re = 100.0
    
    u, v, _, _ = create_test_data(nx, ny)
    
    try:
        u_tp = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
        v_tp = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)
        print("Small grid test: PASS")
    except Exception as e:
        print(f"Small grid test: FAIL - {e}")
        return False
    
    # Test with zero velocity fields
    u_zero = np.zeros_like(u)
    v_zero = np.zeros_like(v)
    
    try:
        u_tp = navier_stokes_x_vectorized(u_zero, v_zero, dx, dy, dt, Re)
        v_tp = navier_stokes_y_vectorized(u_zero, v_zero, dx, dy, dt, Re)
        print("Zero velocity test: PASS")
    except Exception as e:
        print(f"Zero velocity test: FAIL - {e}")
        return False
    
    # Test with very small time step
    dt_small = 1e-10
    
    try:
        u_tp = navier_stokes_x_vectorized(u, v, dx, dy, dt_small, Re)
        v_tp = navier_stokes_y_vectorized(u, v, dx, dy, dt_small, Re)
        print("Small time step test: PASS")
    except Exception as e:
        print(f"Small time step test: FAIL - {e}")
        return False
    
    return True

def test_output_shapes():
    """Test that output arrays have correct shapes."""
    print("\nTesting output shapes...")
    
    nx, ny = 10, 8
    u, v, _, _ = create_test_data(nx, ny)
    dx, dy = 0.1, 0.1
    dt = 0.001
    Re = 100.0
    
    u_tp = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
    v_tp = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)
    
    expected_u_shape = u.shape
    expected_v_shape = v.shape
    
    u_shape_correct = u_tp.shape == expected_u_shape
    v_shape_correct = v_tp.shape == expected_v_shape
    
    print(f"u output shape: {u_tp.shape} (expected: {expected_u_shape}) - {'PASS' if u_shape_correct else 'FAIL'}")
    print(f"v output shape: {v_tp.shape} (expected: {expected_v_shape}) - {'PASS' if v_shape_correct else 'FAIL'}")
    
    return u_shape_correct and v_shape_correct

def main():
    """Run all tests."""
    print("=" * 50)
    print("Testing Vectorized Navier-Stokes Implementation")
    print("=" * 50)
    
    tests = [
        ("Consistency", test_consistency),
        ("Performance", test_performance),
        ("Edge Cases", test_edge_cases),
        ("Output Shapes", test_output_shapes)
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"{test_name} test failed with exception: {e}")
            results.append((test_name, False))
    
    print("\n" + "=" * 50)
    print("Test Summary")
    print("=" * 50)
    
    all_passed = True
    for test_name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"{test_name}: {status}")
        if not result:
            all_passed = False
    
    print(f"\nOverall: {'ALL TESTS PASSED' if all_passed else 'SOME TESTS FAILED'}")
    
    return all_passed

if __name__ == "__main__":
    main() 