class InitialConditions():

    def __init__(self, Re, t, dt, tf, SOR, n_iter):
        self.Re = Re
        self.tf = tf
        self.SOR = SOR
        self.n_iter = n_iter
        self.dt = dt
        self.t = t
        self.N = tf/dt
        self.t = t
        self.it = 0
        self.error = []
