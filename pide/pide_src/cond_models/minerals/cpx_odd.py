#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Zhao2016_WetDiopside(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	a_zhao = [1.29,3.69,3.52]
	r_zhao = [0.58,1.28,0.79]
	e_zhao = [83000.0,115000.0,129000.0]

	if (mechanism == 'dry') or (mechanism == 'polaron') or (mechanism == 'ionic'):
		raise ValueError('Dry (Polaron and Polaron) are not included in electrical conductivity model: Zhao2016_WetDiopside')
	else:
		cond = (a_zhao[0] * (water**(r_zhao[0])) * np.exp(- (e_zhao[0]) / (R_const * T))) +\
		(a_zhao[1] * (water)**(r_zhao[1]) * np.exp(- (e_zhao[1]) / (R_const * T))) +\
		(a_zhao[2] * (water)**(r_zhao[2]) * np.exp(- (e_zhao[2]) / (R_const * T)))

	return cond
	
def Liu2021_WetClinopyroxene_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	#param1 is xFE = feo / (feo+mgo)
	#water has to be converted to wt%

	E = 106000.0
	sigma = 3.44
	n = -0.05
	alpha = 64000.0
	beta = -26000.0
	
	if (mechanism == 'dry') or (mechanism == 'polaron') or (mechanism == 'ionic'):
		raise ValueError('Dry (Polaron and Polaron) are not included in electrical conductivity model: Liu2021_WetClinopyroxene_xFe')
	else:
		cond = (10**sigma) * (xFe**n) * (water) * np.exp(-(E - (alpha * (xFe**(1.0/3.0))) - (beta * (water**(1.0/3.0)))) / (R_const * T))
	
	return cond
	
def Fullea2011_DryCpx_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma_pol_fullea = 10**3.25

	b_ful = [2.075,-2.77,2.61,-1.09] #polynomical coefficients that calculates the fe-dependency of act. enthalpies from Seifert et al (1982)
	fe_pol_fullea = (b_ful[0] + (b_ful[1] * xFe) + (b_ful[2] * (xFe**2.0)) + (b_ful[3]* (xFe**3.0))) * 1e5

	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Fullea2011_DryCpx_xFe')
	else:
		cond = (sigma_pol_fullea * np.exp(-fe_pol_fullea / (R_const * T)))

	return cond