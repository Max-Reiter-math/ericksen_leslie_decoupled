#!/usr/bin/env python
"""
class for a standard benchmark setting for the ericksen-leslie model: 
    a 2D-example where a smooth solution for the ericksen-leslie model exists
"""
from fenics import *
import numpy as np

class smooth_2d:
    # as experiment one but in 2D
    def __init__(self, n=2**5, dt=0.0005, T=2.0, dim=2):
        self.name="smooth 2D"
        if dim != 2: raise ValueError("Dimension "+str(dim)+" is not supported for experiment "+self.name)
        self.dim = dim
        self.fine = n
        self.dt = dt
        self.t0 = 0.0
        self.T = T
        # - model parameters namely v_el, const_A, mu_1, mu_4, mu_5, mu_6, lam
        self.param_dict= {"v_el":1.0, "const_A":1.0, "nu":0.1,"mu_1":1.0, "mu_4": 0.1, "mu_5":1.0, "mu_6":1.0 , "lam":1.0}
        self.parameters = list(self.param_dict.values())
        self.mesh = RectangleMesh(Point(-1,-1),Point(1,1),n,n)
        self.param_dict["dim"]=self.dim
        def boundary(x):
            return x[0] < (DOLFIN_EPS -1) or x[0] > (1 - DOLFIN_EPS) or x[1] < (DOLFIN_EPS -1) or x[1] > (1 - DOLFIN_EPS)
        self.boundary = boundary

        # - initial conditions
        zero_expr = Expression(("0.0","0.0"), degree=2)
        d0_expr = Expression(("sin( 2.0*pi*(cos(x[0])-sin(x[1]) ) )","cos( 2.0*pi*(cos(x[0])-sin(x[1]) ) )"), degree=2, pi = np.pi)
        self.ics = [zero_expr,Constant(0.0),d0_expr,zero_expr]
        # boundary conditions
        self.bcs = [zero_expr, Constant(0.0), d0_expr, zero_expr]

