#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Sun2017_Mudstone_DS6_0_5GPa(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 75000.0
	E2 = 59000.0
	sigma1 = 10.0**1.45
	sigma2 = 10.0**0.92

	if method == 'array':

		T = T[0]
		P = P[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 773.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))

	elif method == 'index':

		if T <= 773.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))


	return cond

def Sun2017_Mudstone_DS7_1_5GPa(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 70000.0
	E2 = 49000.0
	sigma1 = 10.0**1.80
	sigma2 = 10.0**1.18

	if method == 'array':

		T = T[0]
		P = P[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 773.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))

	elif method == 'index':

		if T <= 773.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))


	return cond


def Sun2017_Mudstone_DS8_2_5GPa(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 64000.0
	E2 = 36000.0
	sigma1 = 10.0**1.40
	sigma2 = 10.0**0.94

	if method == 'array':

		T = T[0]
		P = P[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 723.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))

	elif method == 'index':

		if T <= 723.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))


	return cond
	