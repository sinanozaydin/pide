#!/usr/bin/env python3

import numpy as np
from ..eos.fluid_eos import Sanchez_Valle_2013_WaterDensity

R_const = 8.3144621

def Sakuma_2016_NaCl_H2O(T,salinity, P, method = None):

	P = P * 1e3 #converting GPa to MPa

	gamma_11 = [-2.76823e-12, 2.86668e-11, -1.0112e-11]
	gamma_12 = [6.32515e-9, -6.3595e-9, 2.14326e-8]
	gamma_13 = [-2.92588e-6, 2.69121e-5, -9.2074e-6]
	
	gamma_21 = [6.52051e-9, -7.43514e-8, 2.23618e-8]
	gamma_22 = [-1.47966e-5, 1.67038e-4, -4.54299e-5]
	gamma_23 = [6.88977e-3, -7.25629e-2, 1.89836e-2]
	
	gamma_31 = [-2.60077e-6, 3.64027e-5, -7.50611e-6]
	gamma_32 = [6.12874e-3, -9.01143e-2, 1.51621e-2]
	gamma_33 = [-3.17282, 50.2186, -6.22277]
	
	beta_11 = (gamma_11[0] * salinity**2) + (gamma_11[1] * salinity) + gamma_11[2]
	beta_12 = (gamma_12[0] * salinity**2) + (gamma_12[1] * salinity) + gamma_12[2]
	beta_13 = (gamma_13[0] * salinity**2) + (gamma_13[1] * salinity) + gamma_13[2]
	
	beta_21 = (gamma_21[0] * salinity**2) + (gamma_21[1] * salinity) + gamma_21[2]
	beta_22 = (gamma_22[0] * salinity**2) + (gamma_22[1] * salinity) + gamma_22[2]
	beta_23 = (gamma_23[0] * salinity**2) + (gamma_23[1] * salinity) + gamma_23[2]
	
	beta_31 = (gamma_31[0] * salinity**2) + (gamma_31[1] * salinity) + gamma_31[2]
	beta_32 = (gamma_32[0] * salinity**2) + (gamma_32[1] * salinity) + gamma_32[2]
	beta_33 = (gamma_33[0] * salinity**2) + (gamma_33[1] * salinity) + gamma_33[2]
	
	alpha_1 = (beta_11 * T**2) + (beta_12 * T) + beta_13
	alpha_2 = (beta_21 * T**2) + (beta_22 * T) + beta_23
	alpha_3 = (beta_31 * T**2) + (beta_32 * T) + beta_33
	
	cond = (alpha_1 * P**2) + (alpha_2 * P) + alpha_3
	
	return cond

def Bannard_1975_NaCl_H2O(T,salinity,P = None, method = None):

	phi_11 = -1.61994e-12
	phi_12 = 4.32808e-11
	phi_13 = 1.15235e-10
	phi_14 = 2.52257e-10
	
	phi_21 = 1.88235e-9
	phi_22 = -5.82409e-8
	phi_23 = -3.37538e-7
	phi_24 = -4.53779e-7
	
	phi_31 = -5.65158e-7
	phi_32 = 2.70538e-5
	phi_33 = 2.4027e-4
	phi_34 = 2.97574e-4
	
	phi_41 = 4.64690e-5
	phi_42 = -6.70560e-3
	phi_43 = -2.69091e-2
	phi_44 = -8.37212e-2
	
	phi_51 = 2.58834e-3
	phi_52 = 6.92510e-1
	phi_53 = -3.2293
	phi_54 = 9.48091
	
	delta_1 = (phi_11*salinity**3) + (phi_12*salinity**2) + (phi_13*salinity) + phi_14
	delta_2 = (phi_21*salinity**3) + (phi_22*salinity**2) + (phi_23*salinity) + phi_24
	delta_3 = (phi_31*salinity**3) + (phi_32*salinity**2) + (phi_33*salinity) + phi_34
	delta_4 = (phi_41*salinity**3) + (phi_42*salinity**2) + (phi_43*salinity) + phi_44
	delta_5 = (phi_51*salinity**3) + (phi_52*salinity**2) + (phi_53*salinity) + phi_54
	
	cond = (delta_1 * T**4) + (delta_2 * T**3) + (delta_3 * T**2) + (delta_4 * T) + delta_5
	
	return cond
	
def Sinmyo2017(T, P, salinity, method):

	if isinstance(T, (int,float)):
	
		if T > 373.15:
			raise ValueError('The value cannot be lower than 100C(373.15 K) to use Sinmyo2016 model.'+\
				'It might be more appropriate to use another fluid conductivity model for current temperatures.'+\
				'Low temperature model suggestion would be Sakuma_2015_NaCl_H2O')
		else:
			if len(T[T>373.15]) > 0:
				raise ValueError('The values in temperature array cannot be lower than 100C(373.15 K) to use Sinmyo2016 model.' +\
				'It might be more appropriate to use another fluid conductivity model for current temperatures.'+\
				'Low temperature model suggestion would be Sakuma_2015_NaCl_H2O')	
	
	#first calculating pure water density at T and P using iapws08
	
	
	
	if method == 'index':
		rho_water = Sanchez_Valle_2013_WaterDensity(T = T, P = P) / 1e3 #converting to g/cm^3
	else:
		
		P = P[0]
		T = T[0]
		salinity = salinity[0]
		P[P<=0.0] = 1e-4 #Anything below 0 is converted to 1 atm = 1e-4 GPa
		rho_water = Sanchez_Valle_2013_WaterDensity(T = T, P = P) / 1e3
		
	#setting up calculation parameters

	lambda_1 = 1573.0
	lambda_2 = -1212.0
	lambda_3 = 537062.0
	lambda_4 = -208122721.0
	
	def calculate_lambda_0(rho,temp):
		
		_lambda_0 = lambda_1 + (lambda_2 * rho) + (lambda_3 / temp) + (lambda_4 / temp**2)
		
		return _lambda_0
	
	
	lambda_0 = calculate_lambda_0(rho = rho_water, temp = T)
		

	A = -1.7060
	B = -93.78
	C = 0.8075
	D = 3.0781
	
	if method == 'index':
		if salinity > 0.0:
			cond = 10**(A + (B/T) + (C * np.log10(salinity)) + (D * (np.log10(rho_water))) + np.log10(lambda_0))
		elif salinity == 0.0:
			cond = 10**(A + (B/T) + (D * (np.log10(rho_water))) + np.log10(lambda_0))
		else:
			raise ValueError('The salinity value cannot be less than 0.0')
				
	elif method == 'array':
	
		cond = np.zeros(len(T))
		
		for i in range(0,len(T)):
			if salinity[i] > 0.0:
				cond[i] =  10**(A + (B/T[i]) + (C * np.log10(salinity[i])) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))
			elif salinity[i] == 0.0:
				cond[i] = 10**(A + (B/T[i]) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))
			else:
				raise ValueError('The salinity value cannot be less than 0.0')
				

	return cond

def Guo2019(T, P, salinity, method):


	if isinstance(T, (int,float)):
	
		if T > 373.15:
			raise ValueError('The value cannot be lower than 100C(373.15 K) to use Sinmyo2016 model.'+\
				'It might be more appropriate to use another fluid conductivity model for current temperatures.')
		else:
			if len(T[T>373.15]) > 0:
				raise ValueError('The values in temperature array cannot be lower than 100C(373.15 K) to use Sinmyo2016 model.' +\
				'It might be more appropriate to use another fluid conductivity model for current temperatures.')	
				
	#first calculating pure water density at T and P using iapws08

	if method == 'array':

		P = P[0]
		T = T[0]
		salinity = salinity[0]

	rho_water = np.zeros(len(P))

	P[P<=0.0] = 1e-4 #Anything below 0 is converted to 1 atm = 1e-4 GPa
	rho_water = Sanchez_Valle_2013_WaterDensity(T = T, P = P) / 1e3
	
	#setting up calculation parameters

	lambda_1 = 1573.0
	lambda_2 = -1212.0
	lambda_3 = 537062.0
	lambda_4 = -208122721.0

	lambda_0 = lambda_1 + (lambda_2 * rho_water) + (lambda_3 / T) + (lambda_4 / T**2)

	A = -0.919
	B = -872.5
	C = 0.852
	D = 7.61

	if method == 'array':
		
		cond = np.zeros(len(T))
	
		for i in range(0,len(T)):
			if salinity[i] > 0:
				cond[i] = 10**(A + (B/T[i]) + (C * np.log10(salinity[i])) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))
			else:
				cond[i] = 10**(A + (B/T[i]) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))
	else:
		if salinity > 0:
			cond = 10**(A + (B/T) + (C * np.log10(salinity) + (D * (np.log10(rho_water))) + np.log10(lambda_0)))
		else:
			cond = 10**(A + (B/T) + (D * (np.log10(rho_water))) + np.log10(lambda_0))
		
	
	return cond

def Manthilake2021_Aqueous(T, P, salinity, method):

	if method == 'array':
		T = T[0]
	
		for i in range(0,len(T)):
			if T[i] > 673.0:
				cond = (10.0**2.4) * np.exp((-70000.0) / (R_const * T[i]))
			else:
				cond = (10.0**-2.3) * np.exp((-80000.0) / (R_const * T[i]))
	
	else:
		
		if T > 673.0:
			cond = (10.0**2.4) * np.exp((-70000.0) / (R_const * T))
		else:
			cond = (10.0**-2.3) * np.exp((-80000.0) / (R_const * T))
		
	return cond
	

		


