#!/usr/bin/env python
"""
A fixpoint solver for the ericksen-leslie model which makes use of a linear formulation 
where all three variables are decoupled.
"""
from fenics import *
import numpy as np
from .basic_models.basemodel_linear_fp_decoupled import *
from .basic_models.assist_funcs import *

class decoupled_fp_solver(basemodel_linear_fp_decoupled):
    def __init__(self,  submodel,fp_tol = 1e-6, max_iter = 50, silent=True):
        super().__init__(fp_tol=fp_tol, max_iter=max_iter, silent=silent)
        self.submodel = submodel
        self.modelname = "decoupled_fp_solver:"+self.submodel
        self.description = "Decoupled Fixpoint scheme for the general Ericksen-Leslie-model"
        
    def iterate(self):
        """
        in order to list the possible methods for linear solvers use
        list_linear_solver_methods()
        list_krylov_solver_methods()
        list_krylov_solver_preconditioners()
        """
        [bc_v, bc_p, bc_d, bc_q] = self.bcs

        # setting Krylov solver and preconditioner
        Ksolver = KrylovSolver("bicgstab","hypre_amg")
        Ksolver.parameters["absolute_tolerance"] = 1E-10
        Ksolver.parameters["relative_tolerance"] = 1E-9
        Ksolver.parameters["maximum_iterations"] = 1000
        #Ksolver.parameters["linear_solver"] = "mumps"

        # - solving navier-stokes like equation
        Aa = assemble(self.La)
        ba = assemble(self.Ra)
        bc_v.apply(Aa, ba)
        Ksolver.solve(Aa, self.ul.vector(), ba)
        
        # solving the director equation
        Ac = assemble(self.Lc)
        bc = assemble(self.Rc)
        bc_d.apply(Ac, bc)
        solve(Ac, self.dl.vector(), bc, "mumps")
        
        # solving the variational derivative, here the discrete laplacian (q = - \Delta d)
        Ab = assemble(self.Lb)
        bb = assemble(self.Rb)
        bc_q.apply(Ab,bb)
        solve(Ab,self.ql.vector(),bb, "mumps") 
    
    def update_fp_err(self):
        """
        the fixpoint tolerance is measured in terms of the L^2 norm wrt to the velocity and H^1 norm wrt to the director
        """
        (vl,pl)=self.ul.split(deepcopy=True)
        (vl0,pl0)=self.ul0.split(deepcopy=True)

        e_q = self.ql.vector()[:]-self.ql0.vector()[:]
        e_q = np.max(np.abs(e_q))
        e_v = np.sqrt(assemble((vl-vl0)**2*dx))
        e_d = np.sqrt(assemble((self.dl-self.dl0)**2*dx + grad(self.dl-self.dl0)**2*dx))
       
        self.fp_err_list = [e_v, e_d ,e_q]
        self.fp_err = e_v + e_d

    def get_model_dependent_log(self):
        """
        computation of energies
        """
        [v_el, const_A, nu, mu_1, mu_4, mu_5, mu_6, lam] = self.parameters
        (vl,pl)=self.ul.split(deepcopy=True)
        kinetic_energy = 0.5*assemble((vl**2)*dx )
        elastic_energy =  const_A*0.5*assemble(((grad(self.dl))**2)*dx)
        self.energy = [ kinetic_energy + elastic_energy, kinetic_energy, elastic_energy]
        self.energy_labels = ["Total Energy", "Kinetic Energy", "Elastic Energy"]
        return {"energy": dict(zip(self.energy_labels, self.energy))  }
       

    def create_variational_formulation(self):
        """
        choose one of the two submodels
        """
        if self.submodel == "general": 
            self.create_variational_formulation_general()
        elif self.submodel == "simple":
            self.create_variational_formulation_simple()

    def create_variational_formulation_general(self):
        dt = Constant(self.dt)
        [v_el, const_A, nu, mu_1, mu_4, mu_5, mu_6, lam] = self.parameters
        
        # Define variational functions
        (self.vl1,self.pl1) = TrialFunctions(self.TH)
        self.ul0 = Function(self.TH)
        (self.vl0,self.pl0)=split(self.ul0)
        self.ul = Function(self.TH)
        (self.vl,self.pl)=split(self.ul)
        self.u_test = TestFunction(self.TH)
        (self.a,self.h)=split(self.u_test)
        self.u0 = Function(self.TH)
        (self.v0,self.p0)=split(self.u0)
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

        self.init_functions = [self.u0.sub(0),self.u0.sub(1),self.d0,self.q0]

        # projection 
        self.grad_d0_project = Function(self.TensorF) 
        
        # Mass Lumping
        dxL = dx(scheme='vertex', degree=1, metadata={'representation': 'quadrature', 'degree': 1})

        # Energy term
        E = const_A * inner( grad(Constant(0.5) *(self.dl+self.d0)), grad(self.b))

        # Matrix I-outer(d,d)
        I_dd0 = cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim)

        # Leslie stress tensor
        T_L = dt*v_el*Constant(mu_1+lam**2)*inner(inner(0.5 *(self.dl0+self.d0),dot(grad_sym(self.vl1),0.5 *(self.dl0+self.d0))),inner(0.5 *(self.dl0+self.d0),dot(grad_sym(self.a),0.5 *(self.dl0+self.d0))))*dx\
                + dt*Constant(mu_4)*inner( grad_sym(self.vl1), grad_sym(self.a))*dx \
                + v_el*dt* Constant(mu_5+mu_6-lam**2)*inner( dot(grad_sym(self.vl1),self.dl0), dot(grad_sym(self.a),self.dl0))*dx \
                - v_el*dt*Constant(lam)*inner(dot( I_dd0 , self.ql0), dot(grad_sym(self.a),0.5 *(self.dl0+self.d0)))*dxL \
                - v_el*dt*inner(dot(grad_skw(self.a),self.ql0),0.5 *(self.dl0+self.d0))*dxL 

        # Momentum equation (of navier-stokes type) - includes divergence freedom
        # this also means that this incompressible equation system is solved in a coupled way as in the navier stokes case
        # for a coupled solver, however decoupled of the director equation and the discrete laplacian
        Fa =  inner( (self.vl1-self.v0) , self.a )*dx \
            + dt*inner(dot(self.v0, nabla_grad(self.vl1)), self.a)*dx + dt*0.5*div(self.v0)*inner(self.vl1, self.a)*dx \
            - dt*v_el* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.a)) ,  self.ql0 )*dxL\
            -dt*inner(self.pl1,div(self.a))*dx + div(self.vl1)*self.h*dx \
            + T_L

        self.La = lhs(Fa)
        self.Ra = rhs(Fa)

        # discrete Laplacian
        Fb = E*dx - inner(self.ql1,self.b)*dxL 

        self.Lb = lhs(Fb)
        self.Rb = rhs(Fb)

        # director equation
        Fc = inner((self.dl1-self.d0), self.c)*dxL \
            + v_el*dt* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.vl0)) ,  self.c )*dxL \
            + dt*inner(dot(cross_mat(0.5 *(self.dl1+self.d0), 0.5 *(self.dl0+self.d0),self.dim), self.ql0), self.c)*dxL \
            + v_el*dt*Constant(lam)*inner(self.c,dot(cross_mat(0.5 *(self.dl1+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(grad_sym(self.vl0),0.5 *(self.dl0+self.d0))))*dxL \
            - v_el*dt*inner(dot(grad_skw(self.vl0),0.5 *(self.dl1+self.d0)),self.c)*dxL 

        self.Lc = lhs(Fc)
        self.Rc = rhs(Fc)
        
    def create_variational_formulation_simple(self):
        dt = Constant(self.dt)
        [v_el, const_A, nu, mu_1, mu_4, mu_5, mu_6, lam] = self.parameters

        # Define variational functions
        (self.vl1,self.pl1) = TrialFunctions(self.TH)
        self.ul0 = Function(self.TH)
        (self.vl0,self.pl0)=split(self.ul0)
        self.ul = Function(self.TH)
        (self.vl,self.pl)=split(self.ul)
        self.u_test = TestFunction(self.TH)
        (self.a,self.h)=split(self.u_test)
        self.u0 = Function(self.TH)
        (self.v0,self.p0)=split(self.u0)
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

        self.init_functions = [self.u0.sub(0),self.u0.sub(1),self.d0,self.q0]

        #projection
        self.grad_d0_project = Function(self.TensorF)
        
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
            - dt*v_el* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.a)) ,  self.ql0 )*dxL\
            -dt*inner(self.pl1,div(self.a))*dx + div(self.vl1)*self.h*dx 

        self.La = lhs(Fa)
        self.Ra = rhs(Fa)

        # discrete Laplacian
        Fb = E*dx - inner(self.ql1,self.b)*dxL 
        self.Lb = lhs(Fb)
        self.Rb = rhs(Fb)

        # director equation
        Fc = inner((self.dl1-self.d0), self.c)*dxL \
            + dt*v_el* inner( dot(cross_mat(0.5 *(self.dl0+self.d0), 0.5 *(self.dl0+self.d0),self.dim), dot(self.grad_d0_project,self.vl0)) ,  self.c )*dxL \
            + dt*inner(dot(cross_mat(0.5 *(self.dl1+self.d0), 0.5 *(self.dl0+self.d0),self.dim), self.ql0), self.c)*dxL 

        self.Lc = lhs(Fc)
        self.Rc = rhs(Fc)
        