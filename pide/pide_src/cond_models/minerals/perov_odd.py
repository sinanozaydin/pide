#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621
K_const = 8.6173

def Han2025_Perovskite_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	alpha = 3.6
	sigma = 4.8
	E = 0.41
	
	cond = (10**sigma) * (xFe**alpha) * np.exp(-(E - (1e-2 * 0.26 * P)) / (K_const * 1e-5 * T))
	
	return cond
	
	