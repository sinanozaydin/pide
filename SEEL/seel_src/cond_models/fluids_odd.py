#!/usr/bin/env python3

import numpy as np
import iapws

R_const = 8.3144621

def Sinmyo2016(T, P, salinity, method):

	#first calculating pure water density at T and P using iapws08

	P = P * 1e3 #converting GPa to MPa	

	if method == 'index':
		d = iapws.iapws08.SeaWater(T = T, P = P, S = 0)
		rho = d.rho / 1e3 #in g/cm3
		rho_water = rho
	else:
		
		P = P[0]
		T = T[0]
		salinity = salinity[0]
		rho_water = np.zeros(len(P))

		for i in range(0,len(P)):
	
			d = iapws.iapws08.SeaWater(T = T[i], P = P[i], S = 0)
			rho = d.rho / 1e3 #in g/cm3
			rho_water[i] = rho

	#setting up calculation parameters

	lambda_1 = 1573.0
	lambda_2 = -1212.0
	lambda_3 = 537062.0
	lambda_4 = -208122721.0

	lambda_0 = lambda_1 + (lambda_2 * rho) + (lambda_3 / T) + (lambda_4 / T**2)

	A = -1.7060
	B = -93.78
	C = 0.8075
	D = 3.0781

	if method == 'index':
		if salinity > 0:
			cond = 10**(A + (B/T) + (C * np.log10(salinity)) + (D * (np.log10(rho_water))) + np.log10(lambda_0))
		else:
			cond = 10**(A + (B/T) + (D * (np.log10(rho_water))) + np.log10(lambda_0))
	elif method == 'array':
		cond = np.zeros(len(T))
		for i in range(0,len(T)):
			if salinity[i] > 0:
				cond[i] =  10**(A + (B/T[i]) + (C * np.log10(salinity[i])) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))
			else:
				cond[i] = 10**(A + (B/T[i]) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))
		

	return cond

def Guo2019(T, P, salinity, method):

	#first calculating pure water density at T and P using iapws08

	P = P * 1e3 #converting GPa to MPa

	if method == 'array':

		P = P[0]
		T = T[0]
		salinity = salinity[0]

	rho_water = np.zeros(len(P))

	for i in range(0,len(P)):

		d = iapws.iapws08.SeaWater(T = T[i], P = P[i], S = 0)
		rho = d.rho / 1e3 #in g/cm3
		rho_water[i] = rho

	#setting up calculation parameters

	lambda_1 = 1573.0
	lambda_2 = -1212.0
	lambda_3 = 537062.0
	lambda_4 = -208122721.0

	lambda_0 = lambda_1 + (lambda_2 * rho) + (lambda_3 / T) + (lambda_4 / T**2)

	A = -0.919
	B = -872.5
	C = 0.852
	D = 7.61

	for i in range(0,len(T)):
		if salinity[i] > 0:
			cond = 10**(A + (B/T[i]) + (C * np.log10(salinity[i])) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))
		else:
			cond = 10**(A + (B/T[i]) + (D * (np.log10(rho_water[i]))) + np.log10(lambda_0[i]))

	return cond

def Manthilake2021_Aqueous(T, P, salinity, method):

	if method == 'array':
		T = T[0]

	for i in range(0,len(T)):
		if T[i] > 673.0:
			cond = (10.0**2.4) * np.exp((-70000.0) / R_const * T[i])
		else:
			cond = (10.0**-2.3) * np.exp((-80000.0) / R_const * T[i])

	return cond
		


