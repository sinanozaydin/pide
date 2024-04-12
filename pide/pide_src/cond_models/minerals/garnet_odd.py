#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621
	
def Dai2009_DryandWetPyropeGarnet(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A_dai = [1036.0,1950.0]
	B = 0.044
	r_dai_gt = 0.63
	E_dai = [128000.0,70000.0]
	V_dai = [2.50e3,-0.57e3]

	cond_wet = (A_dai[1] * (water**r_dai_gt) * np.exp(-((E_dai[1]) + (P*V_dai[1])) / (R_const * T)))
	cond_dry = (A_dai[0] * (1-(B*P)) * np.exp(-((E_dai[0]) + (P*V_dai[0])) / (R_const * T)))

	if (mechanism == None):
		cond = cond_dry + cond_wet
	elif (mechanism == 'polaron'):
		cond = cond_dry
	elif (mechanism == 'dry'):
		cond = cond_dry
	elif (mechanism == 'proton'):
		cond = cond_wet
	else:
		raise ValueError('Ionic conduction is not included in electrical conductivity model: Dai2009_wetPyropeGarnet')

	return cond
	
def Fullea2011_DryGarnet_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma_pol_fullea = (10.0)**(-0.72 + np.log10(1 - (0.44*P)))
	sigma_i_fullea = 4.96

	e_i_fullea = 205000.0

	b_ful = [2.6,-15.33,80.4,-194.6,202.6,-75.0]
	dv_ful = 2.5 * 1e-6
	
	fe_pol_fullea = (b_ful[0] + (b_ful[1] * xFe) + (b_ful[2] * (xFe**2.0)) + (b_ful[3]* (xFe**3.0)) +\
	(b_ful[4]* (xFe**4.0)) + (b_ful[5]* (xFe**5.0)) + (dv_ful * P)) * 1e5

	cond = (sigma_i_fullea * np.exp(-e_i_fullea / (R_const * T))) +\
	 (sigma_pol_fullea * np.exp(-fe_pol_fullea / (R_const * T)))
	
	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Fullea2011_DryGarnet_xFe')
	else:
		return cond
	
def Liu2021_WetAlmandineGarnet_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	E = 89000.0
	sigma = 4.79
	n = 4.25
	alpha = -4000.0
	beta = -15000.0
		
	if (mechanism == 'proton') or (mechanism == None):
		cond = (10**sigma) * (xFe**n) * (water) * np.exp(-(E - (alpha * (xFe**(1.0/3.0))) - (beta * (water**(1.0/3.0)))) / (R_const * T))
	else:
		raise ValueError('Dry conduction (polaron and ionic) is not included in electrical conductivity model: Liu2021_WetAlmandineGarnet_xFe')
	
	return cond
	
	
	