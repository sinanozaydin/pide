#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Hu2024_Apatite(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	E = 205000.0
	dv = 9.31e3
	sigma = 6.32
	
	cond = 10**sigma * np.exp(-(E + (dv*P)) / (R_const*T))
	
	return cond
	