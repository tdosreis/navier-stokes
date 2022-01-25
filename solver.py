def poisson(i, j, u_tp, v_tp, dx, dy, dt, P, c, SOR):

    cm = (
        u_tp[i, j] - u_tp[i-1, j] + v_tp[i, j] - v_tp[i, j-1]
        )   # conservation of mass

    pressure = (P[i+1, j] + P[i-1, j] + P[i, j+1] + P[i, j-1])

    P[i, j] = SOR * (c[i, j])*(pressure - (dx/dt)*(cm)) + (1.0 - SOR) * P[i, j]

    return P


def streamlines(i, j, dx, dy, dt, phi, vorticity, SOR):

    PHI = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])

    phi[i, j] = (
            (1.0/4.0) * SOR * (PHI + (dx*dy) * vorticity[i, j]) +
            (1.0 - SOR)*phi[i, j]
            )

    return phi
