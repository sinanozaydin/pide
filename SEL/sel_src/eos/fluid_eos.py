#!/usr/bin/env python3

import numpy as np

def Sanchez_Valle_2013_WaterDensity(T,P):
	
	#Function to calculate water density at given T and P
	#T in Kelvin
	
	#Converting P to Pa from GPa
	P = P * 1e9

	a1 = 1.148187e3
	a2 = -2.540804
	a3 = 2.917138e-5
	
	b1 = 8.507742e-3
	b2 = -2.412079e-8
	
	c1 = 1.811854e-11
	c2 = 9.660446e-2
	
	a = a1 + (a2 * T) + (a3 * T**2)
	b = (b1*(P**0.5)) + (b2 * P)
	c = (c1*T*P) + (c2*T*np.log(P))
	
	rho = a + b + c
	
	return rho