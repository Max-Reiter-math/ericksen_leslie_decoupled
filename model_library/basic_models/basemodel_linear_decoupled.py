#!/usr/bin/env python
"""
class that equips basemodel with a space discretization that decouples the three quantities d,q,v
ics refers to inital conditions
bcs refers to boundary conditions
"""
from .basemodel import *
from fenics import *
import numpy as np
import warnings

class basemodel_linear_decoupled(basemodel):
    def __init__(self, silent=True):
        super().__init__(silent=silent)
        self.modelname = "basemodel_linear_decoupled"
    
    def iterate(self):
        """
        will later be overwritten by superclasses
        """
        warnings.warn("method not defined in class", UserWarning)
        pass
       
    def set_ics(self, ics):
        for i in range(len(self.init_functions)):
            assign(self.init_functions[i], interpolate(ics[i], self.init_spaces[i]))
        assign(self.ul0.sub(0), self.u0.sub(0))
        assign(self.ul0.sub(1), self.u0.sub(1))   
        if not self.silent: print("-- projecting d0...")
        self.grad_d0_project.assign(project(grad(self.d0),self.TensorF, solver_type="petsc"))
        # - for consistency also compute q0 - but necessary?
        if not self.silent: print("-- computing q0...")
        Ab = assemble(self.Lb)
        bb = assemble(self.Rb)
        solve(Ab,self.q0.vector(),bb, "mumps") # solver_parameters={'linear_solver': 'mumps'})
        self.dl0.assign(self.d0)
        self.ql0.assign(self.q0)

    def get_ics(self):
        return self.init_functions

    def update_ics(self):
        # sets the current solution as IC
        
        assign(self.u0, self.ul)
        self.d0.assign(self.dl)
        self.q0.assign(self.ql)
        # projection of gradient of d
        self.grad_d0_project.assign(project(grad(self.d0),self.TensorF, solver_type="petsc"))
       
    def set_bcs(self, bcs):
        self.bcs = []
        for i in range(len(self.boundary_spaces)):
            if bcs[i] != None:
                self.bcs.append(DirichletBC(self.boundary_spaces[i], bcs[i], self.boundary))
            else:
                self.bcs.append("")
    
    def get_bcs(self):
        return self.bcs

    def get_functions(self,dc=False):
        (vl0,pl0)=self.ul0.split(deepcopy=dc)
        return [vl0,pl0,self.dl0,self.ql0]

    def create_function_spaces(self, mesh):
        # -- create Function spaces
        # make Taylor-Hood-space for the velocity
        if self.dim==3:
            basic_fem = tetrahedron
        if self.dim==2: 
            basic_fem = triangle
        
        self.V = VectorElement('P', basic_fem, 2 , dim=self.dim)
        self.P = FiniteElement('P', basic_fem, 1 )
        self.element = MixedElement(self.V,self.P)
        self.TH = FunctionSpace(mesh,self.element)
        # - no MixedElement here since quantities are decoupled
        self.D = VectorFunctionSpace(mesh,'P', 1 , dim=self.dim)
        self.Q = VectorFunctionSpace(mesh, 'P', 1 , dim=self.dim)

        # later needed for the projection of grad (d)
        self.TensorF = TensorFunctionSpace(mesh, "P",1, shape=(self.dim,self.dim))

        self.function_spaces = [self.TH,self.D,self.Q,self.TensorF]
        self.boundary_spaces = [self.TH.sub(0),self.TH.sub(1),self.D,self.Q]
        self.init_spaces = [self.TH.sub(0).collapse(),self.TH.sub(1).collapse(),self.D,self.Q]

        self.normal_vector = FacetNormal(mesh)

    def create_variational_formulation(self):
        warnings.warn("method not defined in class", UserWarning)
        pass
        