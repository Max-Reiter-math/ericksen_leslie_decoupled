#!/usr/bin/env python
"""
file that collects all currently available models for an automatic query of those
"""
from experiments import *
from model_library import *

def get_all_experiments():
    # give all experiments with according kwargs
    all_experiments = [[smooth_2d, dict()], [annihilation_2, {"dim":2}], [velocity_driven_flow_2 ,{"dim":2}], [annihilation_2 , {"dim":3}], [velocity_driven_flow_2, {"dim":3}]]  
    return all_experiments

def get_exp_dict():
    all_experiments = {"smooth_2d": smooth_2d, "annihilation_2": annihilation_2, "velocity_driven_flow_2" : velocity_driven_flow_2}
    return all_experiments

def get_model_dict():
    all_models = { "decoupled_fp_solver": decoupled_fp_solver, "decoupled_fp_solver_chorin": decoupled_fp_solver_chorin,}
    return all_models