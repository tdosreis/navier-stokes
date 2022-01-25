import numpy as np
from navier_stokes import navier_stokes_x, navier_stokes_y
from solver import poisson, streamlines
from boundary_conditions import (
    velocity_bc, pressure_bc, pressure_ghosts, vorticity_bc
    )


def simulate_flow(u, v, Nxc, Nyc, dt, tf, n_iter, SOR, Re, Us, Un, Vw, Ve):

    dx = 1.0/Nxc
    dy = 1.0/Nyc
    x = np.linspace(0 - dx, 1.0 + dx, Nxc+2)
    y = np.linspace(0 - dy, 1.0 + dy, Nyc+2)
    N = tf/dt
    error = []
    t = 0   # initial time
    it = 0   # number of iterations

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
    x2d = np.zeros((Nxc+2, Nyc+2))
    y2d = np.zeros((Nxc+2, Nyc+2))
    pressure_bc(Nxc, Nyc, c)
    P_chk = np.zeros((Nxc+2, Nyc+2))

    for n in range(1, int(N+1)):

        for i in range(0, Nxc+1):

            for j in range(0, Nyc+1):

                velocity_bc(i, j, u, v, Us, Un, Vw, Ve, Nxc, Nyc)

        for i in range(1, Nxc):

            for j in range(1, Nyc+1):

                u_tp = navier_stokes_x(i, j, u, u_tp, v, dx, dy, dt, Re)

        for i in range(1, Nxc+1):

            for j in range(1, Nyc):

                v_tp = navier_stokes_y(i, j, u, v_tp, v, dx, dy, dt, Re)

        for it in range(1, n_iter+1):   # SOR

            for i in range(1, Nxc+1):

                for j in range(1, Nyc+1):

                    P_chk[i, j] = P[i, j]

            for i in range(1, Nxc+1):

                for j in range(1, Nyc+1):

                    P = poisson(i, j, u_tp, v_tp, dx, dy, dt, P, c, SOR)

                    phi = streamlines(i, j, dx, dy, dt, phi, vorticity, SOR)

                    pressure_ghosts(Nxc, Nyc, P)

        rmse = np.sqrt(np.sum(np.sum((P_chk - P)**2))/(Nxc*Nyc))

        error.append(rmse)

        for i in range(1, Nxc+1):

            for j in range(1, Nyc+1):

                vorticity_bc(i, j, Nxc, Nyc, dx, dy, phi, vorticity)

        for i in range(1, Nxc):

            for j in range(1, Nyc+1):

                u[i, j] = u_tp[i, j] - (dt/dx)*(P[i+1, j] - P[i, j])

        for i in range(1, Nxc+1):

            for j in range(1, Nyc):

                v[i, j] = v_tp[i, j] - (dt/dy)*(P[i, j+1] - P[i, j])

        x[0:Nxc+2] = np.linspace(0, Nxc, Nxc+2)

        y[0:Nyc+2] = np.linspace(0, Nyc, Nyc+2)

        for i in range(0, Nxc+2):

            for j in range(0, Nyc+2):

                x2d[i, j] = x[i]

                y2d[i, j] = y[j]

        for i in range(0, Nxc+1):

            for j in range(0, Nyc+1):

                u_contour[i, j] = (1.0/2.0) * (u[i, j] + u[i, j+1])

                v_contour[i, j] = (1.0/2.0) * (v[i, j] + v[i+1, j])

                p_contour[i, j] = (
                    (1.0/4.0) * (P[i, j] + P[i+1, j] + P[i, j+1] + P[i+1, j+1])
                    )

                vorticity[i, j] = (
                    (v[i+1, j] - v[i, j])/dx - (u[i, j+1] - u[i, j])/dy
                    )

                phi_contour[i, j] = (
                    (1.0/4.0) * (
                        (phi[i, j] + phi[i+1, j] + phi[i, j+1] + phi[i+1, j+1])
                        )
                    )

        [X, Y] = np.meshgrid(x, y)

        STREAMLINE = np.transpose(phi_contour[0:Nxc+2, 0:Nyc+2])

        t += dt

        it += 1

    return X, Y, STREAMLINE, u_contour, v_contour, p_contour, vorticity
