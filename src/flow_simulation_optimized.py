import numpy as np
from navier_stokes_vectorized import navier_stokes_x_vectorized, navier_stokes_y_vectorized
from solver import poisson, streamlines
from boundary_conditions import (
    velocity_bc, pressure_bc, pressure_ghosts, vorticity_bc
    )


def simulate_flow_optimized(u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve):
    """
    Highly optimized version of flow simulation with maximum vectorization.
    """
    dx = 1.0/Nxc
    dy = 1.0/Nyc
    x = np.linspace(0 - dx, 1.0 + dx, Nxc+2)
    y = np.linspace(0 - dy, 1.0 + dy, Nyc+2)
    N = tf/dt
    error = []
    t = 0
    it = 0

    # Initialize arrays
    c = np.zeros((Nxc+2, Nyc+2))
    u = np.zeros((Nxc+1, Nyc+2))
    u_tp = np.zeros((Nxc+1, Nyc+2))
    v = np.zeros((Nxc+2, Nyc+1))
    v_tp = np.zeros((Nxc+2, Nyc+1))
    P = np.zeros((Nxc+2, Nyc+2))
    u_contour = np.zeros((Nxc+1, Nyc+1))
    v_contour = np.zeros((Nxc+1, Nyc+1))
    p_contour = np.zeros((Nxc+1, Nyc+1))
    phi_contour = np.zeros((Nxc+2, Nyc+2))
    phi = np.zeros((Nxc+2, Nyc+2))
    vorticity = np.zeros((Nxc+2, Nyc+2))
    P_chk = np.zeros((Nxc+2, Nyc+2))
    
    # Set up pressure boundary conditions
    pressure_bc(Nxc, Nyc, c)
    
    # Pre-compute coordinate arrays (only once)
    X, Y = np.meshgrid(x, y, indexing='ij')
    x2d, y2d = X, Y

    for n in range(1, int(N+1)):
        print(f'Running iteration number {n} of {int(N)}', end="\r")

        # OPTIMIZED: Vectorized boundary conditions
        # Set boundary velocities directly
        u[0, :] = Vw  # left boundary
        u[-1, :] = Ve  # right boundary
        v[:, 0] = Us  # bottom boundary  
        v[:, -1] = Un  # top boundary

        # VECTORIZED: Use vectorized Navier-Stokes functions
        u_tp = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
        v_tp = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)

        # SOR iterations for pressure and streamlines
        for it in range(1, n_iter+1):
            # OPTIMIZED: Vectorized array copy
            P_chk[1:Nxc+1, 1:Nyc+1] = P[1:Nxc+1, 1:Nyc+1]

            # Still need loops for SOR (these are inherently sequential)
            for i in range(1, Nxc+1):
                for j in range(1, Nyc+1):
                    P = poisson(i, j, u_tp, v_tp, dx, dy, dt, P, c, SOR)
                    phi = streamlines(i, j, dx, dy, dt, phi, vorticity, SOR)
                    pressure_ghosts(Nxc, Nyc, P)

        # OPTIMIZED: Vectorized RMSE calculation
        rmse = np.sqrt(np.mean((P_chk[1:Nxc+1, 1:Nyc+1] - P[1:Nxc+1, 1:Nyc+1])**2))
        error.append(rmse)

        # OPTIMIZED: Vectorized vorticity boundary conditions
        # Apply boundary conditions for vorticity
        vorticity[0, :] = 0  # left boundary
        vorticity[-1, :] = 0  # right boundary
        vorticity[:, 0] = 0  # bottom boundary
        vorticity[:, -1] = 0  # top boundary

        # OPTIMIZED: Vectorized velocity updates with pressure correction
        u[1:Nxc, 1:Nyc+1] = u_tp[1:Nxc, 1:Nyc+1] - (dt/dx)*(P[2:Nxc+1, 1:Nyc+1] - P[1:Nxc, 1:Nyc+1])
        v[1:Nxc+1, 1:Nyc] = v_tp[1:Nxc+1, 1:Nyc] - (dt/dy)*(P[1:Nxc+1, 2:Nyc+1] - P[1:Nxc+1, 1:Nyc])

        # OPTIMIZED: Vectorized contour calculations
        # Velocity contours
        u_contour[0:Nxc+1, 0:Nyc+1] = 0.5 * (u[0:Nxc+1, 0:Nyc+1] + u[0:Nxc+1, 1:Nyc+2])
        v_contour[0:Nxc+1, 0:Nyc+1] = 0.5 * (v[0:Nxc+1, 0:Nyc+1] + v[1:Nxc+2, 0:Nyc+1])
        
        # Pressure contours
        p_contour[0:Nxc+1, 0:Nyc+1] = 0.25 * (
            P[0:Nxc+1, 0:Nyc+1] + P[1:Nxc+2, 0:Nyc+1] + 
            P[0:Nxc+1, 1:Nyc+2] + P[1:Nxc+2, 1:Nyc+2]
        )
        
        # Vorticity
        vorticity[0:Nxc+1, 0:Nyc+1] = (
            (v[1:Nxc+2, 0:Nyc+1] - v[0:Nxc+1, 0:Nyc+1])/dx - 
            (u[0:Nxc+1, 1:Nyc+2] - u[0:Nxc+1, 0:Nyc+1])/dy
        )
        
        # Stream function contours
        phi_contour[0:Nxc+1, 0:Nyc+1] = 0.25 * (
            phi[0:Nxc+1, 0:Nyc+1] + phi[1:Nxc+2, 0:Nyc+1] + 
            phi[0:Nxc+1, 1:Nyc+2] + phi[1:Nxc+2, 1:Nyc+2]
        )

        # OPTIMIZED: Use pre-computed meshgrid
        STREAMLINE = np.transpose(phi_contour[0:Nxc+2, 0:Nyc+2])

        t += dt
        it += 1

    return X, Y, STREAMLINE, u_contour, v_contour, p_contour, vorticity


def simulate_flow_ultra_optimized(u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve):
    """
    Ultra-optimized version with additional performance improvements.
    """
    dx = 1.0/Nxc
    dy = 1.0/Nyc
    x = np.linspace(0 - dx, 1.0 + dx, Nxc+2)
    y = np.linspace(0 - dy, 1.0 + dy, Nyc+2)
    N = tf/dt
    error = []
    t = 0
    it = 0

    # Initialize arrays with optimal memory layout
    c = np.zeros((Nxc+2, Nyc+2), dtype=np.float64)
    u = np.zeros((Nxc+1, Nyc+2), dtype=np.float64)
    u_tp = np.zeros((Nxc+1, Nyc+2), dtype=np.float64)
    v = np.zeros((Nxc+2, Nyc+1), dtype=np.float64)
    v_tp = np.zeros((Nxc+2, Nyc+1), dtype=np.float64)
    P = np.zeros((Nxc+2, Nyc+2), dtype=np.float64)
    u_contour = np.zeros((Nxc+1, Nyc+1), dtype=np.float64)
    v_contour = np.zeros((Nxc+1, Nyc+1), dtype=np.float64)
    p_contour = np.zeros((Nxc+1, Nyc+1), dtype=np.float64)
    phi_contour = np.zeros((Nxc+2, Nyc+2), dtype=np.float64)
    phi = np.zeros((Nxc+2, Nyc+2), dtype=np.float64)
    vorticity = np.zeros((Nxc+2, Nyc+2), dtype=np.float64)
    P_chk = np.zeros((Nxc+2, Nyc+2), dtype=np.float64)
    
    # Set up pressure boundary conditions
    pressure_bc(Nxc, Nyc, c)
    
    # Pre-compute coordinate arrays and constants
    X, Y = np.meshgrid(x, y, indexing='ij')
    dt_dx = dt/dx
    dt_dy = dt/dy
    dx2 = dx**2
    dy2 = dy**2

    # Pre-allocate slices for better performance
    u_interior = slice(1, Nxc)
    v_interior = slice(1, Nxc+1)
    j_interior_u = slice(1, Nyc+1)
    j_interior_v = slice(1, Nyc)
    
    # Boundary indices
    u_left, u_right = 0, -1
    v_bottom, v_top = 0, -1

    for n in range(1, int(N+1)):
        print(f'Running iteration number {n} of {int(N)}', end="\r")

        # OPTIMIZED: Direct boundary assignment
        u[u_left, :] = Vw
        u[u_right, :] = Ve
        v[:, v_bottom] = Us
        v[:, v_top] = Un

        # VECTORIZED: Use vectorized Navier-Stokes functions
        u_tp = navier_stokes_x_vectorized(u, v, dx, dy, dt, Re)
        v_tp = navier_stokes_y_vectorized(u, v, dx, dy, dt, Re)

        # SOR iterations (still sequential due to dependencies)
        for it in range(1, n_iter+1):
            P_chk[1:Nxc+1, 1:Nyc+1] = P[1:Nxc+1, 1:Nyc+1]

            for i in range(1, Nxc+1):
                for j in range(1, Nyc+1):
                    P = poisson(i, j, u_tp, v_tp, dx, dy, dt, P, c, SOR)
                    phi = streamlines(i, j, dx, dy, dt, phi, vorticity, SOR)
                    pressure_ghosts(Nxc, Nyc, P)

        # OPTIMIZED: Vectorized error calculation
        error.append(np.sqrt(np.mean((P_chk[1:Nxc+1, 1:Nyc+1] - P[1:Nxc+1, 1:Nyc+1])**2)))

        # OPTIMIZED: Vectorized boundary conditions and updates
        # Set vorticity boundaries to zero
        vorticity[0, :] = vorticity[-1, :] = vorticity[:, 0] = vorticity[:, -1] = 0

        # OPTIMIZED: Vectorized velocity updates with pre-computed constants
        u[u_interior, j_interior_u] = (
            u_tp[u_interior, j_interior_u] - 
            dt_dx * (P[2:Nxc+1, j_interior_u] - P[1:Nxc, j_interior_u])
        )
        
        v[v_interior, j_interior_v] = (
            v_tp[v_interior, j_interior_v] - 
            dt_dy * (P[v_interior, 2:Nyc+1] - P[v_interior, 1:Nyc])
        )

        # OPTIMIZED: Vectorized contour calculations with pre-computed slices
        u_contour[:] = 0.5 * (u[0:Nxc+1, 0:Nyc+1] + u[0:Nxc+1, 1:Nyc+2])
        v_contour[:] = 0.5 * (v[0:Nxc+1, 0:Nyc+1] + v[1:Nxc+2, 0:Nyc+1])
        
        p_contour[:] = 0.25 * (
            P[0:Nxc+1, 0:Nyc+1] + P[1:Nxc+2, 0:Nyc+1] + 
            P[0:Nxc+1, 1:Nyc+2] + P[1:Nxc+2, 1:Nyc+2]
        )
        
        vorticity[0:Nxc+1, 0:Nyc+1] = (
            (v[1:Nxc+2, 0:Nyc+1] - v[0:Nxc+1, 0:Nyc+1])/dx - 
            (u[0:Nxc+1, 1:Nyc+2] - u[0:Nxc+1, 0:Nyc+1])/dy
        )
        
        phi_contour[0:Nxc+1, 0:Nyc+1] = 0.25 * (
            phi[0:Nxc+1, 0:Nyc+1] + phi[1:Nxc+2, 0:Nyc+1] + 
            phi[0:Nxc+1, 1:Nyc+2] + phi[1:Nxc+2, 1:Nyc+2]
        )

        STREAMLINE = np.transpose(phi_contour[0:Nxc+2, 0:Nyc+2])
        t += dt
        it += 1

    return X, Y, STREAMLINE, u_contour, v_contour, p_contour, vorticity 