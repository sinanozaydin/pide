#!/usr/bin/env python3

import numpy as np
R_const = 8.3144621
	
def Rauch2002_OpxSol(T,P,depth,h2o_fug,o2_fug,fe_opx,al_opx,method):

	h2o_fug = h2o_fug * 1e5 #converting to barsfac
	P = P * 1e9 #converting to GPa to Pa

	A_rauch = 0.01354 #ppm/bar
	E_rauch = -4563.0 #J/mol
	dv = 12.1e-6 #converted to m3/mol
	
	max_opx_water = A_rauch * h2o_fug * np.exp(-(E_rauch) / (R_const*T)) * np.exp(-(dv * P) / (R_const*T)) 
	
	return max_opx_water
	
def Mierdel2007_OpxSol(T,P,depth,h2o_fug,o2_fug,fe_opx,al_opx,method):

	water_non_al = Rauch2002_OpxSol(T = T, P = P, depth = depth,h2o_fug=h2o_fug,o2_fug = o2_fug, fe_opx = fe_opx, al_opx = al_opx,method = method)[0]
	
	A = 0.042 #ppm/bar
	h2o_fug = h2o_fug * 1e5 #converting to bar
	E = -79.685e3 #J/mol
	dv = 11.3e-6 #m3/mol
	P = P * 1e9 #converting to GPa
	
	water_al = A * (h2o_fug**0.5) * np.exp(-E / (R_const*T)) * np.exp(-(dv * P) / (R_const*T))
	
	return water_non_al + water_al
	
	
def Kang2022_OpxSol(T,P,depth,h2o_fug,o2_fug,fe_opx,al_opx,method):

	water_non_al = Rauch2002_OpxSol(T = T, P = P, depth = depth,h2o_fug=h2o_fug,o2_fug = o2_fug, fe_opx = fe_opx, al_opx = al_opx,method = method)[0]
	
	#water fugacity is in GPa
	
	B2 = 1800 #ppm/GPa
	E2 = 9.5e3 #j/mol
	dv2 = 8.2e-6 #m3/mol
	P = P * 1e9 #converting to Pa for pressure dependence
	
	max_opx_water = water_non_al + (B2 * al_opx * (h2o_fug**(0.5)) * np.exp(-(E2 + (dv2 * P)) / (R_const*T)))
	
	return max_opx_water
	
def Guo2020_OpxSol(T,P,depth,h2o_fug,o2_fug,fe_opx,al_opx,method):

	A = 0.00263 #ppm/bar
	dv = 17.04e-6 #m3/mol
	P = P * 1e9 #converting to Pa or J/m3
	h2o_fug = h2o_fug * 1e5 #converting GPa to bar
	n = 1.24
	
	max_opx_water = A * (h2o_fug**n) * np.exp(-(P * dv) / (R_const*T))
	
	return max_opx_water
