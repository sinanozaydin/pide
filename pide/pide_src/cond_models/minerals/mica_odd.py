#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Li2016_Phlogopite(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A_li_001 = 10**10.15
	A_li_010 = 10**8.41
	A_li_110 = 10**6.09

	E_li_001 = 204000.0
	E_li_010 = 179000.0
	E_li_110 = 134000.0

	cond_001 = A_li_001 * np.exp(-(E_li_001) / (R_const * T))
	cond_010 = A_li_010 * np.exp(-(E_li_010) / (R_const * T))
	cond_110 = A_li_110 * np.exp(-(E_li_110) / (R_const * T))

	cond = (cond_001 * cond_010 * cond_110)**(1.0/3.0)

	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity models of micas')
	else:
		return cond

def Li2017_Phlogopite_param1_F(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	#param1 is F content in weight percent
	
	if np.any(param1 == 0.0):
		raise ValueError('param1 for mica has to be set to a value for this function. It denotes to F content in phlogopite.')
	
	A_li = 10**8.59
	E_li = 191000.0
	r_li = 0.98
	
	cond = A_li * (param1**r_li) * np.exp(-E_li / (R_const * T))

	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity models of micas')
	else:
		return cond
