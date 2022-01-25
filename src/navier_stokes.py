def navier_stokes_x(i, j, u, u_tp, v, dx, dy, dt, Re):

    uu_f = ((1.0/2.0) * (u[i+1, j] + u[i, j]))**2.0
    uu_b = ((1.0/2.0) * (u[i-1, j] + u[i, j]))**2.0

    Duu_x = (uu_f - uu_b)/dx

    u_sum_f = u[i, j] + u[i, j+1]
    u_sum_b = u[i, j] + u[i, j-1]

    v_sum_f = v[i, j] + v[i+1, j]
    v_sum_b = v[i, j-1] + v[i+1, j-1]

    uv_x_f = (1.0/2.0)*(u_sum_f)*(1.0/2.0)*(v_sum_f)
    uv_x_b = (1.0/2.0)*(u_sum_b)*(1.0/2.0)*(v_sum_b)

    D2u_x2 = ((u[i+1, j] - 2.0*u[i, j] + u[i-1, j])/dx**2.0)
    D2u_y2 = ((u[i, j+1] - 2.0*u[i, j] + u[i, j-1])/dy**2.0)

    Duv_y = (uv_x_f - uv_x_b)/dy

    advection_u = Duu_x + Duv_y
    diffusion_u = D2u_x2 + D2u_y2

    u_tp[i, j] = u[i, j] + dt * (-advection_u + (1.0/Re) * diffusion_u)

    return u_tp


def navier_stokes_y(i, j, u, v_tp, v, dx, dy, dt, Re):

    vv_f = ((1.0/2.0) * (v[i, j+1] + v[i, j]))**2.0
    vv_b = ((1.0/2.0) * (v[i, j] + v[i, j-1]))**2.0

    Dvv_y = (vv_f - vv_b)/dy

    u_sum_f = u[i, j+1] + u[i, j]
    u_sum_b = u[i-1, j+1] + u[i-1, j]

    v_sum_f = v[i, j] + v[i+1, j]
    v_sum_b = v[i, j] + v[i-1, j]

    uv_y_f = (1.0/2.0)*(u_sum_f)*(1.0/2.0)*(v_sum_f)
    uv_y_b = (1.0/2.0)*(u_sum_b)*(1.0/2.0)*(v_sum_b)

    D2v_x2 = ((v[i+1, j] - 2.0*v[i, j] + v[i-1, j])/dx**2.0)
    D2v_y2 = ((v[i, j+1] - 2.0*v[i, j] + v[i, j-1])/dy**2.0)

    Duv_x = (uv_y_f - uv_y_b)/dx

    advection_v = Duv_x + Dvv_y
    diffusion_v = D2v_x2 + D2v_y2

    v_tp[i, j] = v[i, j] + dt * (-advection_v + (1.0/Re) * diffusion_v)

    return v_tp
