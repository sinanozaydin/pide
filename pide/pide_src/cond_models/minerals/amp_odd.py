#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Shen2021(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	E = 52210.0
	dv = 0.33e3
	sigma = 10.0**1.93
	
	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Shen2021')
	else:
		cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond
	
def Shen2020_TremoliteDehydrogenated(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	E = 39210.0
	dv = 0.41e3
	sigma = 10.0**0.73
	
	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Shen2020_TremoliteDehydrogenated')
	else:
		cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond
