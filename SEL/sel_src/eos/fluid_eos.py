#!/usr/bin/env python3

import numpy as np
import math
from scipy import optimize

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
	
def Pitzer_and_Sterner_1994_PureWaterEOS(T, P):

	#Water fugacity calculation for pure water from the EOS of Pitzer and Sterner (1994), adapted after the python script of Tony Withers taken from his personal website.

	coeff = [[0.0,0.0,0.24657688e6,0.51359951e2,0.0,0.0],[0.0,0.0,0.58638965e0,-0.28646939e-2,0.31375577e-4,0.0],
	[0.0,0.0,-0.62783840e1,0.14791599e-1,0.35779579e-3,0.15432925e-7],[0.0,0.0,0,-0.42719875e0,-0.16325155e-4,0.0],
	[0.0,0.0,0.56654978e4,-0.16580167e2,0.76560762e-1,0],[0,0,0,0.10917883e0,0,0.0],
	[0.38878656e13,-0.13494878e9,0.30916564e6,0.75591105e1,0,0],[0,0,-0.65537898e5,0.18810675e3,0.0,0.0],
	[-0.14182435e14,0.18165390e9,-0.19769068e6,-0.23530318e2,0.0,0.0],[0.0,0.0,0.92093375e5,0.12246777e3,0.0,0.0]]

	c = []

	def PSeos(volume, temperature, targetP):  # cc/mol, Kelvins, bars
		R=8314510  # Pa.cc/K/mol
		den=1/volume  # mol/cc

		for i in range(10):
				c.insert(i,coeff[i][0]*temperature**-4+coeff[i][1]*temperature**-2
						+coeff[i][2]*temperature**-1+coeff[i][3]
						+coeff[i][4]*temperature+coeff[i][5]*temperature**2)

		pressure = (den+c[0]*den**2-den**2*((c[2]+2*c[3]*den+3*c[4]*den**2
				+4*c[5]*den**3)/(c[1]+c[2]*den+c[3]*den**2+c[4]*den**3
				+c[5]*den**4)**2)+c[6]*den**2*math.exp(-c[7]*den)
				+c[8]*den**2*math.exp(-c[9]*den))*R*temperature/1e5
		return pressure-targetP  # bars

	def PSvolume(pressure, temperature):  # bars, Kelvins

		volume = optimize.root(PSeos, 10, args = (temperature, pressure))
		return volume.x

	def PSfugacity(pressure, temperature):  # bars, Kelvins

		for i in range(10):
				c.insert(i,coeff[i][0]*temperature**-4+coeff[i][1]*temperature**-2
						+coeff[i][2]*temperature**-1+coeff[i][3]
						+coeff[i][4]*temperature+coeff[i][5]*temperature**2)

		volume=PSvolume(pressure, temperature)
		R = 8314510  # Pa.cc/K/mol
		den = 1/volume  # mol/cc
		fug = math.exp(math.log(den)+c[0]*den+(1/(c[1]+c[2]*den+c[3]*den**2
					+c[4]*den**3+c[5]*den**4)-1/c[1])
					-c[6]/c[7]*(math.exp(-c[7]*den)-1)
					-c[8]/c[9]*(math.exp(-c[9]*den)-1)
					+pressure*1e5/(den*R*temperature)
					+math.log(R*temperature)-1)/1e5
		return fug  # bars

	h2o_fug = np.zeros(len(T))

	for i in range(0,len(T)):

		h2o_fug[i] = PSfugacity(P[i]*1e4,float(T[i])) / 1e4 #in GPa
			
	return h2o_fug
