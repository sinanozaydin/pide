#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621
K_const = 8.6173

def Han2025_Perovskite_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	alpha = 3.81
	sigma = 5.44
	E = 90000.0
	beta = 92000.0 
	
	cond = (10**sigma) * (xFe**alpha) * np.exp(-(E - (beta * (xFe**0.333333))) / (K_const * T))
	
	return cond
	
	