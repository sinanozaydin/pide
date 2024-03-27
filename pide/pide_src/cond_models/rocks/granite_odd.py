#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Guo2017_MuscoviteGranite(T,P,water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 92000.0
	E2 = 216000.0
	sigma1 = 10.0**1.0
	sigma2 = 10.0**7.88

	if method == 'array':

		T = T[0]
		P = P[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 919.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))

	elif method == 'index':

		if T <= 919.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))

	
	return cond

def Guo2017_BiotiteGranite(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 48000.0
	E2 = 206000.0
	sigma1 = 10.0**-1.47
	sigma2 = 10.0**6.68

	if method == 'array':

		T = T[0]
		P = P[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 990.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))

	elif method == 'index':

		if T <= 919.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))

	
	return cond
