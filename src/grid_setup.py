import numpy as np
from initial_conditions import InitialConditions
from boundary_conditions import BoundaryConditions


class SetupGrid(BoundaryConditions, InitialConditions):

    def __init__(self, Re, SOR, n_iter, t, dt, tf, Un, Us, Vw, Ve, Nxc, Nyc):
        InitialConditions.__init__(self, Re, t, dt, tf, SOR, n_iter)
        BoundaryConditions.__init__(self, Nxc, Nyc, Un, Us, Vw, Ve)
        self.c = np.zeros((Nxc+2, Nyc+2))
        self.u = np.zeros((Nxc+1, Nyc+2))
        self.u_tp = np.zeros((Nxc+1, Nyc+2))
        self.v = np.zeros((Nxc+2, Nyc+1))
        self.v_tp = np.zeros((Nxc+2, Nyc+1))
        self.P = np.zeros((Nxc+2, Nyc+2))
        self.u_contour = np.zeros((Nxc+1, Nyc+1))
        self.v_contour = np.zeros((Nxc+1, Nyc+1))
        self.p_contour = np.zeros((Nxc+1, Nyc+1))
        self.phi_contour = np.zeros((Nxc+2, Nyc+2))
        self.phi = np.zeros((Nxc+2, Nyc+2))
        self.vorticity = np.zeros((Nxc+2, Nyc+2))
        self.x2d = np.zeros((Nxc+2, Nyc+2))
        self.y2d = np.zeros((Nxc+2, Nyc+2))
        self.P_chk = np.zeros((Nxc+2, Nyc+2))
