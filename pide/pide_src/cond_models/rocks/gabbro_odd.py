#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2015_DryGabbro(T,P,water, param1, fo2 = None, fo2_ref = None, method = None):

	#the behaviour of sigma is pretty linear, so I interpolated the p-sigma as what to put on the graph.
	sigma_through_p = np.array([1.76,1.67,1.52,1.38])
	p_study = np.array([0.5,1,1.5,2])

	sigma_interp = np.interp(P, p_study, sigma_through_p) #interpolation

	sigma = 10.0**sigma_interp
	E = 72e3
	dv = -6.8e3 #to be multiplied by GPa

	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond

def Wang2022_DryGabbro_param1_xCpx(T, P, water, param1 = None, fo2 = None, fo2_ref = None, method = None):
	
	if param1 is None:
		
		raise ValueError('xCpx has to be defined between 0.1 to 0.9 to calculate use Wang2022_DryGabbro_param1_xCPx. Define param1 as xCpx.')
		
	C = 10**2.47
	D = 125e3
	beta = -8.31*1e3
	gamma = 1.6
	alpha = 0.12
	
	xCpx = param1
	
	cond = C * (xCpx**alpha) * np.exp(-(D + (beta * (xCpx**gamma)))/ (R_const*T))
	
	return cond

	
