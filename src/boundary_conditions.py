import numpy as np

class BoundaryConditions(): 
    
    def __init__(self, Nxc, Nyc, Un, Us, Vw, Ve): 
        self.Nxc = Nxc
        self.Nyc = Nyc
        self.Un = Un
        self.Us = Us
        self.Vw = Vw
        self.Ve = Ve
        self.dx = 1.0/Nxc
        self.dy = 1.0/Nyc
        self.x = np.linspace(0 - self.dx, 1.0 + self.dx, Nxc+2)
        self.y = np.linspace(0 - self.dy, 1.0 + self.dy, Nyc+2)
        
    def pressure_bc(self, c):
        c[2:self.Nxc, 2:self.Nyc] = 1.0/4.0
        c[1, 1] = 1.0/2.0
        c[self.Nxc, 1] = 1.0/2.0
        c[1, self.Nyc] = 1.0/2.0
        c[self.Nxc, self.Nyc] = 1.0/2.0
        c[2:self.Nxc, 1] = 1.0/3.0
        c[1, 2:self.Nyc] = 1.0/3.0
        c[self.Nxc, 2:self.Nyc] = 1.0/3.0
        c[2:self.Nxc, self.Nyc] = 1.0/3.0
        return c

    def velocity_bc(self, i, j, u, v):
        u[i, 0] = 2.0 * self.Us - u[i, 1]
        u[i, self.Nyc+1] = 2.0 * self.Un - u[i, self.Nyc]
        u[0, j+1] = 0.0
        u[self.Nxc, j+1] = u[0, j+1]
        v[0, j] = 2.0 * self.Vw - v[1, j]
        v[self.Nxc+1, j] = 2.0 * self.Ve - v[self.Nxc, j]
        v[i+1, 0] = 0.0
        v[i+1, self.Nyc] = v[i+1, 0]
        return v, u

    def pressure_ghosts(self, P):
        P[0:self.Nxc+2, 0] = 0
        P[0:self.Nxc+2, self.Nyc+1] = 0
        P[0, 0:self.Nyc+2] = 0
        P[self.Nxc+1, 0:self.Nyc+2] = 0
        return P

    def vorticity_bc(self, i, j, phi, vorticity):
        vorticity[i, 0] = - 2.0 * phi[i, 1]/(self.dx*self.dy)
        vorticity[i, self.Nyc + 1] = -2.0 * phi[i, self.Nyc]/(self.dx * self.dy) - 2.0/self.dx
        vorticity[0, j] = -2.0 * phi[1, j]/(self.dx*self.dy)
        vorticity[self.Nxc+1, j] = - 2.0 * phi[self.Nxc, j]/(self.dx*self.dy)
        return vorticity
