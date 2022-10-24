#!/usr/bin/env python
"""
class for a standard benchmark setting for the ericksen-leslie model: 
    annihilation of two defects without an initial flow 
"""
from fenics import *
import numpy as np

class annihilation_2:
    def __init__(self, dim = 3, n=2**4, dt=0.0005, T= 0.1):
        self.name="annihilation"
        if dim in [2,3]: 
            self.dim = dim
        else: 
            raise ValueError("Dimension "+str(dim)+" is not supported for experiment "+self.name)
        self.fine = n
        self.dt = dt
        self.t0 = 0
        self.T = T
        # - model parameters namely v_el, const_A, nu, mu_1, mu_4, mu_5, mu_6, lam
        self.param_dict= {"v_el":.25, "const_A":1.0, "nu":1.0,"mu_1":1.0, "mu_4": 1.0, "mu_5":1.0, "mu_6":1.0 , "lam":1.0}
        self.parameters = list(self.param_dict.values())
        self.param_dict["dim"]=self.dim
        if self.dim == 3:
            self.mesh = BoxMesh(Point(-0.5, -0.5,-0.5), Point(0.5, 0.5,0.5),n,n,n)
            self.boundary = boundary_3d
            zero_expr = Expression(("0.0","0.0","0.0"), degree=2)
            d0_expr = d0_expr_3d
            bc_expr = Expression(("4.0*pow(x[0],2)+4*pow(x[1],2)-0.25","2.0*x[1]","0"),degree=1)
        else:
            self.mesh = RectangleMesh(Point(-0.5,-0.5),Point(0.5,0.5),n,n)
            self.boundary = boundary_2d
            zero_expr = Expression(("0.0","0.0"), degree=2)
            d0_expr = d0_expr_2d
            bc_expr = Expression(("4.0*pow(x[0],2)+4*pow(x[1],2)-0.25","2.0*x[1]"),degree=1)
        # - initial conditions
        self.ics = [zero_expr,Constant(0.0),d0_expr(),zero_expr]
        # boundary conditions
        self.bcs = [zero_expr, Constant(0.0), bc_expr,zero_expr]
    
def boundary_3d(x):
        return x[0] < (DOLFIN_EPS -0.5) or x[0] > (0.5 - DOLFIN_EPS) or x[1] < (DOLFIN_EPS -0.5) or x[1] > (0.5 - DOLFIN_EPS) or x[2] < (DOLFIN_EPS -0.5) or x[2] > (0.5 - DOLFIN_EPS)  

def boundary_2d(x):
        return x[0] < (DOLFIN_EPS -0.5) or x[0] > (0.5 - DOLFIN_EPS) or x[1] < (DOLFIN_EPS -0.5) or x[1] > (0.5 - DOLFIN_EPS)

class d0_expr_3d(UserExpression):
    # - define custom Expression
    def eval(self,values,x):
        if x[1]==0 and (x[0]==0.25 or x[0]==-0.25):
            values[0]=0.0
            values[1]=0.0
            values[2]=1.0
        else:
            values[0]=4.0*pow(x[0],2)+4*pow(x[1],2)-0.25
            values[1]=2.0*x[1]
            values[2]=0.0
            tmp_abs = np.sqrt(values[0]**2+values[1]**2+values[2]**2)
            values[0]= values[0] / tmp_abs
            values[1]= values[1] / tmp_abs
            values[2]= values[2] / tmp_abs
    def value_shape(self):
        return (3,)

class d0_expr_2d(UserExpression):
    # - define custom Expression
    def eval(self,values,x):
        eta = 0.05
        values[0]=pow(x[0],2)+pow(x[1],2)-0.25
        values[1]=x[1]
        tmp_abs = values[0]**2+values[1]**2
        if tmp_abs > DOLFIN_EPS:
            values[0]= values[0] / np.sqrt(tmp_abs + eta**2)
            values[1]= values[1] / np.sqrt(tmp_abs + eta**2)
    def value_shape(self):
        return (2,)