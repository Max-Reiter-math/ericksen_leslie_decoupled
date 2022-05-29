#!/usr/bin/env python
"""
class that equips basemodel with a fixpoint solver
"""
from .basemodel import *
import warnings
import time

class basemodel_fp(basemodel):
    def __init__(self,  fp_tol = 1e-6, max_iter = 50, silent=True):
        super().__init__( silent)
        self.fp_tol = fp_tol
        self.max_iter = max_iter
        self.modelname = "basemodel_fp"

    def do_time_step(self):
        """
        standard fixpoint iteration which can later be specified by overwriting the error measure self.update_fp_err()
        """
        iteration=0
        self.fp_err = float("inf")
        while self.fp_err > self.fp_tol and iteration<self.max_iter:
            iteration+= 1
            if not self.silent: print("time-step ",self.t, ", iteration ",iteration)
            time_meas_start = time.process_time()
                       
            
            self.iterate()
            self.update_fp_err()
            # important to update the fp_error first, bc of information loss after overwriting with update
            self.update()
            
            
            time_meas_end = time.process_time()
            self.computation_time += time_meas_end - time_meas_start
            if not self.silent: print("fp error: ",self.fp_err)

        if self.fp_err<=self.fp_tol:
            return (iteration, self.fp_err_list)
        elif iteration>=self.max_iter: 
            raise RuntimeError("Fixpoint scheme failed to converge in "+str(iteration)+" iterations..")
            return (iteration, self.fp_err_list)

