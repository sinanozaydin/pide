#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2005_DryEnstatite(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None):

	sigma1 = 10**3.78
	E = 171500.0
	dv = 0.03e3

	cond = sigma1 *  np.exp(-(E + (P*dv)) / (R_const * T))

	return cond
	
def Zhang2016_DryOpx_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None):

	#Conductivity model from Zhang et al. (2016, Contr. Min. & Petr.)

	R_const = 8.3144621

	a_z1 = 855610
	e_z1 = 251000.0

	dv_z1 = 4.15*1e3
	a_z2 = 163
	e_z2 = 233000.0

	alpha_z2 = 199000.0

	dv_z2 = 1.06*1e3
	beta_z2 = 12000

	cond = a_z1 * np.exp(- (e_z1 + (P * dv_z1)) / (R_const * T)) +\
	 (a_z2 * xFe * np.exp(- (e_z2 - (alpha_z2 * (xFe**0.33)) + (P * (dv_z2 - (beta_z2*xFe)))) / (R_const * T)))

	return cond
	
def Fullea2011_DryOpx_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None):

	R_const = 8.3144621

	sigma_pol_fullea = 10**3.0 #average value of used models Fullea et al. (2011)

	b_ful = [1.9,-2.77,2.61,-1.09] #polynomical coefficients that calculates the fe-dependency of act. enthalpies from Seifert et al (1982)
	fe_pol_fullea = (b_ful[0] + (b_ful[1] * xFe) + (b_ful[2] * (xFe**2.0)) + (b_ful[3]* (xFe**3.0))) * 1e5

	cond = (sigma_pol_fullea * np.exp(-fe_pol_fullea / (R_const * T)))

	return cond
	
