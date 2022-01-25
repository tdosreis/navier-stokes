import numpy as np
from navier_stokes import navier_stokes_x, navier_stokes_y
from solver import poisson, streamlines
from grid_setup import SetupGrid
import itertools


class FlowSimulation(SetupGrid):

    def __init__(self, Re, SOR, n_iter, t, dt, tf, Un, Us, Vw, Ve, Nxc, Nyc):
        self.Re = Re
        self.tf = tf
        self.SOR = SOR
        self.n_iter = n_iter
        self.dt = dt
        self.t = t
        self.Un = Un
        self.Us = Us
        self.Vw = Vw
        self.Ve = Ve
        self.Nxc = Nxc
        self.Nyc = Nyc
        super().__init__(Re, SOR, n_iter, t, dt, tf, Un, Us, Vw, Ve, Nxc, Nyc)

    def simulate_flow(self):

        self.c = self.pressure_bc(self.c)

        for n in range(1, int(self.N+1)):
            
            print(f'Running iteration number {n} of {int(self.N)}', end="\r")

            x_ = np.arange(0, self.Nxc+1, 1)
            y_ = np.arange(0, self.Nyc+1, 1)

            for (i, j) in itertools.product(x_, y_):

                self.velocity_bc(i, j,
                                 self.u,
                                 self.v)

            x_ = np.arange(1, self.Nxc, 1)
            y_ = np.arange(1, self.Nyc+1, 1)

            for (i, j) in itertools.product(x_, y_):

                self.u_tp = navier_stokes_x(i, j,
                                            self.u,
                                            self.u_tp,
                                            self.v,
                                            self.dx,
                                            self.dy,
                                            self.dt,
                                            self.Re)

            x_ = np.arange(1, self.Nxc+1, 1)
            y_ = np.arange(1, self.Nyc, 1)

            for (i, j) in itertools.product(x_, y_):

                self.v_tp = navier_stokes_y(i, j,
                                            self.u,
                                            self.v_tp,
                                            self.v,
                                            self.dx,
                                            self.dy,
                                            self.dt,
                                            self.Re)

            for it in range(1, self.n_iter+1):

                x_ = np.arange(1, self.Nxc+1, 1)
                y_ = np.arange(1, self.Nyc+1, 1)

                for (i, j) in itertools.product(x_, y_):

                    self.P_chk[i, j] = self.P[i, j]

                x_ = np.arange(1, self.Nxc+1, 1)
                y_ = np.arange(1, self.Nyc+1, 1)

                for (i, j) in itertools.product(x_, y_):

                    self.P = poisson(i, j,
                                     self.u_tp,
                                     self.v_tp,
                                     self.dx,
                                     self.dy,
                                     self.dt,
                                     self.P,
                                     self.c,
                                     self.SOR)

                    self.phi = streamlines(i, j,
                                           self.dx,
                                           self.dy,
                                           self.dt,
                                           self.phi,
                                           self.vorticity,
                                           self.SOR)

                    self.pressure_ghosts(self.P)

            x_ = np.arange(1, self.Nxc+1, 1)
            y_ = np.arange(1, self.Nyc+1, 1)

            for (i, j) in itertools.product(x_, y_):

                self.vorticity_bc(i, j,
                                  self.phi,
                                  self.vorticity)

            x_ = np.arange(1, self.Nxc, 1)
            y_ = np.arange(1, self.Nyc+1, 1)

            for (i, j) in itertools.product(x_, y_):

                self.u[i, j] = (
                    self.u_tp[i, j] -
                    (self.dt/self.dx)*(self.P[i+1, j] - self.P[i, j])
                )

            x_ = np.arange(1, self.Nxc+1, 1)
            y_ = np.arange(1, self.Nyc, 1)

            for (i, j) in itertools.product(x_, y_):

                self.v[i, j] = (
                    self.v_tp[i, j] -
                    (self.dt/self.dy)*(self.P[i, j+1] - self.P[i, j])
                )

            self.x[0:self.Nxc+2] = np.linspace(0, self.Nxc, self.Nxc+2)
            self.y[0:self.Nyc+2] = np.linspace(0, self.Nyc, self.Nyc+2)

            x_ = np.arange(0, self.Nxc+2, 1)
            y_ = np.arange(0, self.Nyc+2, 1)

            for (i, j) in itertools.product(x_, y_):

                self.x2d[i, j] = self.x[i]
                self.y2d[i, j] = self.y[j]

            x_ = np.arange(0, self.Nxc+1, 1)
            y_ = np.arange(0, self.Nyc+1, 1)

            for (i, j) in itertools.product(x_, y_):

                self.u_contour[i, j] = (
                    (1.0/2.0) * (self.u[i, j] + self.u[i, j+1])
                )

                self.v_contour[i, j] = (
                    (1.0/2.0) * (self.v[i, j] + self.v[i+1, j])
                )

                self.p_contour[i, j] = (
                    (1.0/4.0) * (self.P[i, j] + self.P[i+1, j] +
                                 self.P[i, j+1] + self.P[i+1, j+1])
                    )

                self.vorticity[i, j] = (
                    (self.v[i+1, j] - self.v[i, j])/self.dx -
                    (self.u[i, j+1] - self.u[i, j])/self.dy
                    )

                self.phi_contour[i, j] = (
                    (1.0/4.0) * (self.phi[i, j] + self.phi[i+1, j] +
                                 self.phi[i, j+1] + self.phi[i+1, j+1])
                    )

            [self.X, self.Y] = np.meshgrid(self.x, self.y)

            self.STREAMLINE = (
                np.transpose(
                    self.phi_contour[0:self.Nxc+2, 0:self.Nyc+2]
                )
            )

            self.t += self.dt

            self.it += 1

        return (self.X,
                self.Y,
                self.STREAMLINE,
                self.u_contour,
                self.v_contour,
                self.p_contour,
                self.vorticity)
