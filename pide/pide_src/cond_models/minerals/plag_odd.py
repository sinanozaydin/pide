#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Hu2022_WetPlagioclase_Na_param1(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	r = 1.47
	sigma_wet = 9139.0
	E_wet = 85e3
	
	sigma1 = 13.23
	sigma2 = 7.3
	E1 = 128e3
	E2 = 117e3
	alpha1 = -0.56
	alpha2 = -0.59
	beta1 = -1.83
	beta2 = -0.59
	
	tcrit = 873.0
	
	if param1.any() == 0.0:
		raise ValueError('param1 has to be set as Na content. The Na content values cannot be zero for model Hu2022_WetPlagioclase_Na_param1')
	
	if method == 'array':

		T = T[0]
		P = P[0]
		param1 = param1[0]
		water = water[0]
		
		cond_dry = np.zeros(len(T))
		cond_wet = np.zeros(len(T))

		for i in range(0,len(T)):
			if mechanism != 'proton':
			
				if T[i] <= tcrit:
					cond_dry[i] = (sigma1 * (param1[i]**beta1) *  np.exp(-(E1 + (alpha1 * param1[i])) / (R_const * T[i])))
					
				else:
					cond_dry[i] = (sigma2 * (param1[i]**beta2) *  np.exp(-(E2 + (alpha2 * param1[i])) / (R_const * T[i])))
				
			cond_wet[i] = (sigma_wet * (water[i]**r) * np.exp(-E_wet / (R_const * T[i])))
			
	elif method == 'index':
		
		if mechanism != 'proton':
			if T <= tcrit:
				cond_dry = (sigma1 * (param1**beta1) *  np.exp(-(E1 + (alpha1 * param1)) / (R_const * T)))
			else:
				cond_dry = (sigma2 * (param1**beta2) *  np.exp(-(E2 + (alpha2 * param1)) / (R_const * T)))

		cond_wet = (sigma_wet * (water**r) * np.exp(-E_wet / (R_const * T)))
	
	if mechanism == 'proton':
		
		cond = cond_wet
	
	elif (mechanism == 'polaron') or (mechanism == None):
	
		cond = cond_dry
		
	else:
		
		cond = cond_dry + cond_wet
		
	return cond
	
def Hu2015_DryAnorthite(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma_list = np.array([4.62,4.43,4.36])
	P_list = np.array([1,2,3])
	#Sigma seems to be linear with pressure, so fitting a curve to make pressure dependence.
	
	E = 183e3
	dv = 2.39e3
	
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
		
	if (mechanism == 'proton') or (mechanism == 'ionic'):
	
		raise ValueError('Proton and ionic conduction is not included in electrical conductivity model: Hu2015_DryAnorthite')
		
	else:
			
		return cond
	
def Hu2011_DryAlbite(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma_list = np.array([4,3.92,3.89])
	P_list = np.array([1,2,3])
	#Sigma seems to be linear with pressure, so fitting a curve to make pressure dependence.
	
	E = 82e3
	dv = 1.45e3
	
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
		
	if (mechanism == 'proton') or (mechanism == 'ionic'):
	
		raise ValueError('Proton and ionic conduction is not included in electrical conductivity model: Hu2011_DryAlbite')
		
	else:
			
		return cond
	

	