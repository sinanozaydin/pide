#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Hu2014_DryOrthoclase(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma_list = np.array([4.13,4.06,4.00])
	P_list = np.array([1,2,3])
	#Sigma seems to be linear with pressure, so fitting a curve to make pressure dependence.
	
	E = 98e3
	dv = 1.46e3
	
	if method == 'array':

		T = T[0]
		P = P[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):
			
			sigma_interp = np.interp(P[i], P_list, sigma_list)
			
			cond[i] = (10**sigma_interp) * np.exp(-(E + (P[i] * dv)) / (R_const * T[i]))
			
	elif method == 'index':
	
		sigma_interp = np.interp(P, P_list, sigma_list)
			
		cond = (10**sigma_interp) * np.exp(-(E + (P * dv)) / (R_const * T))
		
	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Hu2014_DryOrthoclase')
	else:
		return cond

	
def Dai2018_DryOrthoclase_001(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	
	sigma_list = np.array([4.678,4.69,4.72])
	P_list = np.array([1,2,3])
	#Sigma seems to be linear with pressure, so fitting a curve to make pressure dependence.
	
	E = 104e3
	dv = 2.51e3
	
	if method == 'array':

		T = T[0]
		P = P[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):
			
			sigma_interp = np.interp(P[i], P_list, sigma_list)
			
			cond[i] = (10**sigma_interp) * np.exp(-(E + (P[i] * dv)) / (R_const * T[i]))
			
	elif method == 'index':
	
		sigma_interp = np.interp(P, P_list, sigma_list)
			
		cond = (10**sigma_interp) * np.exp(-(E + (P * dv)) / (R_const * T))
		
	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Dai2018_DryOrthoclase_001')
	else:
		return cond