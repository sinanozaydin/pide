#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Yoshino2012b_DryWad(T, P, water, xFe ,param1, param2, fo2 = None, fo2_ref = None, method = None):

	A = 5.51 #S/m
	E = 94000.0 #j/mol
	dv = -0.62e-6 #converting go m3/mol
	P = P * 1e9 #converting to Pa / J/m3
	
	cond = A * np.exp(-(E + (P*dv)) / (R_const*T))
	
	return cond