def pressure_bc(Nxc, Nyc, c):

    c[2:Nxc, 2:Nyc] = 1.0/4.0
    c[1, 1] = 1.0/2.0
    c[Nxc, 1] = 1.0/2.0
    c[1, Nyc] = 1.0/2.0
    c[Nxc, Nyc] = 1.0/2.0
    c[2:Nxc, 1] = 1.0/3.0
    c[1, 2:Nyc] = 1.0/3.0
    c[Nxc, 2:Nyc] = 1.0/3.0
    c[2:Nxc, Nyc] = 1.0/3.0

    return c


def velocity_bc(i, j, u, v, Us, Un, Vw, Ve, Nyc, Nxc):

    u[i, 0] = 2.0 * Us - u[i, 1]
    u[i, Nyc+1] = 2.0 * Un - u[i, Nyc]
    u[0, j+1] = 0.0
    u[Nxc, j+1] = u[0, j+1]
    v[0, j] = 2.0 * Vw - v[1, j]
    v[Nxc+1, j] = 2.0 * Ve - v[Nxc, j]
    v[i+1, 0] = 0.0
    v[i+1, Nyc] = v[i+1, 0]

    return v, u


def pressure_ghosts(Nxc, Nyc, P):

    P[0:Nxc+2, 0] = 0
    P[0:Nxc+2, Nyc+1] = 0
    P[0, 0:Nyc+2] = 0
    P[Nxc+1, 0:Nyc+2] = 0

    return P


def vorticity_bc(i, j, Nxc, Nyc, dx, dy, phi, vorticity):

    vorticity[i, 0] = - 2.0 * phi[i, 1]/(dx*dy)
    vorticity[i, Nyc + 1] = -2.0 * phi[i, Nyc]/(dx * dy) - 2.0/dx
    vorticity[0, j] = -2.0 * phi[1, j]/(dx*dy)
    vorticity[Nxc+1, j] = - 2.0 * phi[Nxc, j]/(dx*dy)

    return vorticity
