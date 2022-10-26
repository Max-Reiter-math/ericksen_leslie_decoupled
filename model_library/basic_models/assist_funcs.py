#!/usr/bin/env python
"""
File consists of all definitions which are used for the definition of the variational formulation.
"""
from fenics import *
import numpy as np

def grad_sym(v):
    """
    symmetric part of the gradient
    """
    return   0.5*(grad(v) + grad(v).T)
def grad_skw(v):
    """
    skew-symmetric part of the gradient
    """
    return   0.5*(grad(v) - grad(v).T)

def cross_mat(d0, d1, dim):
    """
    dimension independent decomposition of the matrix I-outer(d,d)
    """
    if dim ==2:
        S = as_matrix( [ [d1[1] , d1[0]] , [-d1[0] , d1[1]]] )
        D = as_matrix( [ [1.0 , 0.0] , [0.0 , 0.0]] )
        Sinv = as_matrix( [ [d0[1] , -d0[0]] , [d0[0] , d0[1]]] )
        return (S*D*Sinv)
    if dim ==3:
        I= as_matrix( [ [1.0 ,0.0 , 0.0],  [0.0 , 1.0, 0.0] , [0.0 , 0.0, 1.0]] )
        return (I - outer(d1,d0))

