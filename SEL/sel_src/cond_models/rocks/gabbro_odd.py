#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2015_DryGabbro(T,P,water, param1, param2, fo2 = None, fo2_ref = None, method = None):

	#the behaviour of sigma is pretty linear, so I interpolated the p-sigma as what to put on the graph.
	sigma_through_p = np.array([1.76,1.67,1.52,1.38])
	p_study = np.array([0.5,1,1.5,2])

	sigma_interp = np.interp(P, p_study, sigma_through_p) #interpolation

	sigma = 10.0**sigma_interp
	E = 72e3
	dv = -6.8e3

	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond

def Wang2022(T, P, water, param1, param2, fo2 = None, fo2_ref = None, method = None):

	sigma = 92.40 * (1 - (0.23 * P))
	E = 102e3
	dv = 6.0e3

	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond

	
