#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Zhao2016_WetDiopside(T, P, water, param1, param2, fo2 = None, fo2_ref = None, method = None):

	h2o = water / (1e4) #converting ppm to wt

	a_zhao = [1.29,3.69,3.52]
	r_zhao = [0.58,1.28,0.79]
	e_zhao = [83000.0,115000.0,129000.0]

	cond = (a_zhao[0] * (h2o**(r_zhao[0])) * np.exp(- (e_zhao[0]) / (R_const * T))) +\
	(a_zhao[1] * (h2o)**(r_zhao[1]) * np.exp(- (e_zhao[1]) / (R_const * T))) +\
	(a_zhao[2] * (h2o)**(r_zhao[2]) * np.exp(- (e_zhao[2]) / (R_const * T)))

	return cond
	
def Liu2021_WetOmphacite_param1_Fe(T, P, water, param1, param2, fo2 = None, fo2_ref = None, method = None):

	#param1 is xFE = feo / (feo+mgo)
	#water has to be converted to wt%

	E = 106000.0
	sigma = 3.44
	n = -0.05
	alpha = 64000.0
	beta = -26000.0
	
	water = water / 1e4 #converting it to wt%
	
	cond = (10**sigma) * (param1**n) * (water) * np.exp(-(E - (alpha * (param1**(1.0/3.0))) - (beta * (water**(1.0/3.0)))) / (R_const * T))
	
	return cond
	
def Fullea2011_DryCpx_param1_Fe(T, P, water, param1, param2, fo2 = None, fo2_ref = None, method = None):

	sigma_pol_fullea = 10**3.25

	b_ful = [2.075,-2.77,2.61,-1.09] #polynomical coefficients that calculates the fe-dependency of act. enthalpies from Seifert et al (1982)
	fe_pol_fullea = (b_ful[0] + (b_ful[1] * param1) + (b_ful[2] * (param1**2.0)) + (b_ful[3]* (param1**3.0))) * 1e5

	cond = (sigma_pol_fullea * np.exp(-fe_pol_fullea / (R_const * T)))

	return cond