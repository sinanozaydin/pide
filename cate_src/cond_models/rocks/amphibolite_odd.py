#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Shen2021(T, P, water, method):

	E = 52210.0
	dv = 0.33e3
	sigma = 10.0**1.93
	
	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond

def Wang2012_PH_Perpendicular(T,P,water,method):

	E1 = 66e3
	E2 = 333e3
	sigma1 = 10.0**1.18
	sigma2 = 10.0**17.99

	if method == 'array':

		T = T[0]
		P = P[0]
		
		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 824.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))
	
	elif method == 'index':

		if T <= 824.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))

	return cond

def Wang2012_PH_Parallel(T,P,water,method):

	E1 = 64e3
	E2 = 319e3
	sigma1 = 10.0**1.49
	sigma2 = 10.0**18.43

	if method == 'array':

		T = T[0]
		P = P[0]
		
		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 771.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))
	
	elif method == 'index':

		if T <= 771.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))

	return cond

def Wang2012_HA(T,P,water,method):

	E1 = 67e3
	E2 = 378e3
	sigma1 = 10.0**2.03
	sigma2 = 10.0**23.88

	if method == 'array':

		T = T[0]
		P = P[0]
		
		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 771.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))
	
	elif method == 'index':

		if T <= 736.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))

	return cond
	