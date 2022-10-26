#!/usr/bin/env python
"""
basemodel that comprises the most basic functions of a FEM model for better readability
and better structuring
"""
class basemodel:
    def __init__(self, silent=True):
        self.silent = silent
        self.computation_time = 0.0
        self.modelname = "basemodel"
    
    def init_experiment(self, experiment):
        self.parameters = experiment.parameters
        self.dim = experiment.dim
        self.dt = experiment.dt
        self.dh = experiment.fine
        self.t = experiment.t0
        self.T = experiment.T
        self.boundary = experiment.boundary
        self.mesh = experiment.mesh
        self.exp_ics = experiment.ics
        self.exp_bcs = experiment.bcs
        self.init_discretization()

    def init_discretization(self):
        if not self.silent: print("- creating function spaces...")
        self.create_function_spaces(self.mesh)
        if not self.silent: print("- creating the variational formulation...")
        self.create_variational_formulation()
        if not self.silent: print("- setting the initial conditions...")
        self.set_ics(self.exp_ics)
        if not self.silent: print("- setting the boundary conditions...")
        self.set_bcs(self.exp_bcs)

    def do_time_step(self):
        self.iterate()
        self.update_ics()