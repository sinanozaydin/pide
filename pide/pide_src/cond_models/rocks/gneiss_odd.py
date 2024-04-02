#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2018_Gneiss_DS12(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 63e3
	E2 = 78e3
	dv1 = -7.1e3
	dv2 = -2.69e3
	sigma1 = 10.0**-0.2
	sigma2 = 10.0**1.11

	if method == 'array':

		T = T[0]
		P = P[0]
		
		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 773.0:
				cond[i] = sigma1 * np.exp(-(E1 + (P[i] * dv1)) / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-(E2 + (P[i] * dv2)) / (R_const * T[i]))
	
	elif method == 'index':

		if T <= 773.0:
			cond = sigma1 * np.exp(-(E1 + (P * dv1)) / (R_const * T))
		else:
			cond = sigma2 * np.exp(-(E2 + (P * dv2)) / (R_const * T))

	return cond

def Dai2018_Gneiss_DS13(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 35e3
	E2 = 84e3
	sigma1 = 10.0**-0.92
	sigma2 = 10.0**2.26

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


def Dai2018_Gneiss_DS14(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	E1 = 38e3
	E2 = 87e3
	sigma1 = 10.0**-0.49
	sigma2 = 10.0**2.63

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

def Yu2011_Gneiss_J85(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):
	
	E1 = 56e3
	E2 = 421e3
	sigma1 = 10.0**-0.15
	sigma2 = 10.0**20.88

	if method == 'array':

		T = T[0]
		P = P[0]
		
		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 881.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))
	
	elif method == 'index':

		if T <= 881.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))

	return cond

def Yu2011_Gneiss_J88(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):
	
	E1 = 31e3
	E2 = 224e3
	sigma1 = 10.0**-0.95
	sigma2 = 10.0**11.37

	if method == 'array':

		T = T[0]
		P = P[0]
		
		cond = np.zeros(len(T))

		for i in range(0,len(T)):

			if T[i] <= 799.0:
				cond[i] = sigma1 * np.exp(-E1 / (R_const * T[i]))
			else:
				cond[i] = sigma2 * np.exp(-E2 / (R_const * T[i]))
	
	elif method == 'index':

		if T <= 799.0:
			cond = sigma1 * np.exp(-E1 / (R_const * T))
		else:
			cond = sigma2 * np.exp(-E2 / (R_const * T))

	return cond
