#!/usr/bin/env python
"""
A fixpoint solver for the ericksen-leslie model which makes use of a linear formulation 
where all three variables are decoupled.
For this model conditional existence and conditional convergence against a dissipative solution is known. 
Evtl. a time-step control is necessary.
TO-DO: 
time-step control
"""
from fenics import *
import numpy as np
from .basic_models.basemodel_linear_fp_decoupled_chorin import *
from .basic_models.assist_funcs import *

class decoupled_fp_solver_chorin(basemodel_linear_fp_decoupled_chorin):
    def __init__(self,  submodel,fp_tol = 1e-6, max_iter = 50, silent=True):
        super().__init__(fp_tol=fp_tol, max_iter=max_iter, silent=silent)
        self.submodel = submodel
        self.modelname = "decoupled_fp_solver_chorin:"+self.submodel
        self.description = "Fully decoupled Fixpoint scheme using a chorin projection for the general Ericksen-Leslie-model"
        
    def iterate(self):
        """
        in order to list the possible methods for linear solvers use
        list_linear_solver_methods()
        list_krylov_solver_methods()
        list_krylov_solver_preconditioners()
        """
        [bc_v, bc_p, bc_d, bc_q] = self.bcs

        # setting Krylov solver and preconditioner
        # well known settings for navier-stokes
        Ksolver = KrylovSolver("bicgstab","hypre_amg")
        Ksolver.parameters["absolute_tolerance"] = 1E-10
        Ksolver.parameters["relative_tolerance"] = 1E-9
        Ksolver.parameters["maximum_iterations"] = 1000
        #Ksolver.parameters["linear_solver"] = "mumps"

        # - solving navier-stokes like equation
        Aa = assemble(self.La)
        ba = assemble(self.Ra)
        bc_v.apply(Aa, ba)
        solve(Aa, self.vl.vector(), ba, "gmres", "default")

        Aa2 = assemble(self.La2)
        ba2 = assemble(self.Ra2)
        bc_p.apply(Aa2, ba2)
        solve(Aa2, self.pl.vector(), ba2)

        Aa3 = assemble(self.La3)
        ba3 = assemble(self.Ra3)
        bc_v.apply(Aa3, ba3)
        solve(Aa3, self.vl.vector(), ba3, "gmres", "default")
        
        # solving the director equation
        Ac = assemble(self.Lc)
        bc = assemble(self.Rc)
        solve(Ac, self.dl.vector(), bc, "mumps")
        
        # solving compatibility condition - equation for the variational derivative q
        Ab = assemble(self.Lb)
        bb = assemble(self.Rb)
        solve(Ab,self.ql.vector(),bb, "mumps") # solver_parameters={'linear_solver': 'mumps'})
    
    def update_fp_err(self):
        """
        using L-infty norm as fixpoint tolerance
        """
        e_v = self.vl.vector()[:]-self.vl0.vector()[:]
        e_d = self.dl.vector()[:]-self.dl0.vector()[:]
        e_q = self.ql.vector()[:]-self.ql0.vector()[:]
        
        self.fp_err_list = [np.max(np.abs(e_v)), np.max(np.abs(e_d)) ,np.max(np.abs(e_q))]
        self.fp_err = np.max(self.fp_err_list)    

    def get_model_dependent_log(self):
        """
        computation of metrics depending on this concrete time scheme
        """
        kinetic_energy = assemble((self.vl**2)/2*dx )
        elastic_energy =  assemble((Constant(self.parameters[1])*(grad(self.dl))**2)/2*dx)
        self.energy = [ kinetic_energy + elastic_energy, kinetic_energy, elastic_energy]
        self.energy_labels = ["Total Energy", "Kinetic Energy", "Elastic Energy"]
        return {"energy": dict(zip(self.energy_labels, self.energy))  }
       
    def update_time_step_size(self):
        """
        updates the variable dt if a time step control is activated
        the coefficient in the fixpoint iteration should only depend on time-independent factors 
        and the gradient of the director of the preceeding iteration
        """
        quotient = assemble( ((grad(self.d0))**2)*dx) / assemble( ((grad(self.dl))**2)*dx)
        self.dt = self.dt * quotient


    def create_variational_formulation(self):
        """
        To-do: move the variable declarations for the two submodels which coincide into this method
        """
        # Define variational functions
        self.vl1 = TrialFunction(self.V)
        self.pl1 = TrialFunction(self.P)
        self.vl = Function(self.V)
        self.pl = Function(self.P)
        self.vl0 = Function(self.V)
        self.pl0 = Function(self.P)
        self.v0 = Function(self.V)
        self.p0 = Function(self.P)
        self.a = TestFunction(self.V)
        self.q = TestFunction(self.P)
        #
        self.c = TestFunction(self.D)
        self.d0 = Function(self.D)
        #
        self.b = TestFunction(self.Q)
        self.q0 = Function(self.Q)
        #
        self.dl1 = TrialFunction(self.D)
        self.dl0 = Function(self.D)
        self.dl = Function(self.D)
        #
        self.ql1 = TrialFunction(self.Q)
        self.ql0 = Function(self.Q)
        
        self.ql = Function(self.Q)

        self.init_functions = [self.v0,self.p0,self.d0,self.q0]


        #projection 
        self.grad_d0_project = Function(self.TensorF) #solver_parameters={'linear_solver': 'mumps'})
        if self.submodel == "general": 
            self.create_variational_formulation_general()
        elif self.submodel == "simple":
            self.create_variational_formulation_simple()

    def create_variational_formulation_general(self):
        dt = Constant(self.dt)
        # the following is fine, although its a copy since the constants dont change
        [v_el, const_A, nu, mu_1, mu_4, mu_5, mu_6, lam] = self.parameters
        
        
        # Mass Lumping
        dxL = dx(scheme='vertex', degree=1, metadata={'representation': 'quadrature', 'degree': 1})

        # Energy term
        E = const_A * inner( grad(Constant(0.5) *(self.dl+self.d0)), grad(self.b))

        # Matrix I-outer(d,d)
        I_dd0 = cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim)

        # Leslie stress tensor
        T_L = dt*v_el*Constant(mu_1+lam**2)*inner(inner(0.5 *(self.dl0+self.d0),dot(grad_sym(self.vl1),0.5 *(self.dl0+self.d0))),inner(0.5 *(self.dl0+self.d0),dot(grad_sym(self.a),0.5 *(self.dl0+self.d0))))*dx\
                + dt*Constant(mu_4)*inner( grad_sym(self.vl1), grad_sym(self.a))*dx \
                + dt*v_el* Constant(mu_5+mu_6-lam**2)*inner( dot(grad_sym(self.vl1),self.dl0), dot(grad_sym(self.a),self.dl0))*dx \
                - dt*Constant(lam)*inner(dot( I_dd0 , self.ql0), dot(grad_sym(self.a),0.5 *(self.dl0+self.d0)))*dxL \
                - dt*inner(dot(grad_skw(self.a),self.ql0),0.5 *(self.dl0+self.d0))*dxL 

        # Momentum equation (of navier-stokes type) - includes divergence freedom
        # this also means that this incompressible equation system is solved in a coupled way as in the navier stokes case
        # for a coupled solver, however decoupled of the director equation and the compatibility condition
        Fa =  inner( (self.vl1-self.v0) , self.a )*dx \
            + dt*inner(dot(self.v0, nabla_grad(self.vl1)), self.a)*dx + dt*0.5*div(self.v0)*inner(self.vl1, self.a)*dx \
            - dt*v_el* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.a)) ,  self.ql0 )*dxL\
            + T_L

        self.La = lhs(Fa)
        self.Ra = rhs(Fa)

        # -- Chorin projection steps
        # - pressure correction
        self.La2 = inner(grad(self.pl1),grad(self.q))*dx 
        self.Ra2 = -(1/dt) * div(self.vl)*self.q*dx

        # - velocity correction
        self.La3 = inner(self.vl1,self.a)*dx
        self.Ra3 = inner(self.vl,self.a)*dx - dt* inner(grad(self.pl),self.a)*dx

        # compatibility condition / equation for the variational derivative
        Fb = E*dx - inner(self.ql1,self.b)*dxL 
        self.Lb = lhs(Fb)
        self.Rb = rhs(Fb)

        # director equation
        Fc = inner((self.dl1-self.d0), self.c)*dxL \
            + dt*v_el* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.vl0)) ,  self.c )*dxL \
            + dt*inner(dot(cross_mat(0.5 *(self.dl1+self.d0), 0.5 *(self.dl0+self.d0),self.dim), self.ql0), self.c)*dxL \
            + dt*v_el*Constant(lam)*inner(self.c,dot(cross_mat(0.5 *(self.dl1+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(grad_sym(self.vl0),0.5 *(self.dl0+self.d0))))*dxL \
            - dt*v_el*inner(dot(grad_skw(self.vl0),0.5 *(self.dl1+self.d0)),self.c)*dxL 

        self.Lc = lhs(Fc)
        self.Rc = rhs(Fc)
        
    def create_variational_formulation_simple(self):
        dt = Constant(self.dt)
        # the following is fine, although its a copy since the constants dont change
        [v_el, const_A, nu, mu_1, mu_4, mu_5, mu_6, lam] = self.parameters

        

        #projection
        self.grad_d0_project = Function(self.TensorF) #solver_parameters={'linear_solver': 'mumps'})
        
        # Mass Lumping
        dxL = dx(scheme='vertex', degree=1, metadata={'representation': 'quadrature', 'degree': 1})

        # Energy term
        E = const_A * inner( grad(Constant(0.5) *(self.dl+self.d0)), grad(self.b))

        # Momentum equation (of navier-stokes type) - includes divergence freedom
        # this also means that this incompressible equation system is solved in a coupled way as in the navier stokes case
        # for a coupled solver, however decoupled of the director equation and the compatibility condition
        Fa =  inner( (self.vl1-self.v0) , self.a )*dx \
            + Constant(nu)*dt*inner(grad(self.vl1), grad(self.a))*dx \
            + dt*inner(dot(self.v0, nabla_grad(self.vl1)), self.a)*dx + dt*0.5*div(self.v0)*inner(self.vl1, self.a)*dx \
            - dt*v_el* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.a)) ,  self.ql0 )*dxL

        self.La = lhs(Fa)
        self.Ra = rhs(Fa)
        
        # -- Chorin projection steps
        # - pressure correction
        self.La2 = inner(grad(self.pl1),grad(self.q))*dx 
        self.Ra2 = -(1/dt) * div(self.vl)*self.q*dx

        # - velocity correction
        self.La3 = inner(self.vl1,self.a)*dx
        self.Ra3 = inner(self.vl,self.a)*dx - dt* inner(grad(self.pl),self.a)*dx

        # compatibility condition / equation for the variational derivative
        Fb = E*dx - inner(self.ql1,self.b)*dxL 
        self.Lb = lhs(Fb)
        self.Rb = rhs(Fb)

        # director equation
        Fc = inner((self.dl1-self.d0), self.c)*dxL \
            + dt*v_el* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.vl0)) ,  self.c )*dxL \
            + dt*inner(dot(cross_mat(0.5 *(self.dl1+self.d0), 0.5 *(self.dl0+self.d0),self.dim), self.ql0), self.c)*dxL 

        self.Lc = lhs(Fc)
        self.Rc = rhs(Fc)
        