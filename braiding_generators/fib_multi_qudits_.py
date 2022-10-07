#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 08-31-2022

@author: Abdellah Tounsi
******

'fib_multi_qudits' module computes elementary braiding generators of Fibonacci
(SU(2)_3) anyonic model.

Elementary braiding generators (sigma_n) are computed for any number of anyons 
grouped in qudits. The general form of the state considered in this module is
illustrated by the following example:
    Example:
        1 1 1 1 1 1 1 1 1
        \/  / \/  / \/  /
        i\ /  k\ /  e\ /
          \     /     /
          j\  l/     /f
            \ /     / 
            m\     /
              \   /
               \ /
               t|
        =[(|((1, 1)_i, 1)_j| (X) |((1, 1)_k, 1)_l|)_m (X) |((1, 1)_e, 1)_f|]_t
        state is represented by Python dict 
        {'qudits': [[i, j], [k, l], [e, f]], 'roots': [m, t]}

TODO:
    - raise ValueError's
    - Translate to Cpp
"""

import pickle
from copy import deepcopy
from cplot import cplot

def find_basis(n_qudits, qudit_len):
    """
    generates all states that form the basis of Hilbert space
    of anyons grouped by qudits and fused qudit by qudit.
    
    Inputs:
        n_qudits: int:
            number of qudits.
        qudit_len: int:
            number of outcomes representing one qudit.
    """
    
    f = open(f"./braiding_generators/n{n_qudits}a{qudit_len}/basis.pkl","rb")
    states = pickle.load(f)
    f.close()
    
    return states
 
def braiding_generator(index, n_qudits, qudit_len, show=True):
    """
    calculates matrix representation of the braiding generator -in the basis
    of multi-qudit fusion space- which exchanges
    index'th anyon with the (index + 1)'th anyon.
    
    Inputs:
        index: int:
            index of braiding operator.
        n_qudits: int:
            number of qudits.
        qudit_len: int:
            number of outcomes representing one qudit.
    Returns:
        (numpy.array whose dimension equals to the dimension of
        anyons' Hilbert space, basis)
    """
    f = open(f"./braiding_generators/n{n_qudits}a{qudit_len}/braiding_generator{index}.pkl","rb")
    
    sig = pickle.load(f)
    f.close()
    
    return sig