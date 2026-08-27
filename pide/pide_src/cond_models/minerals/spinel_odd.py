#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2009_DryPyropeGarnet_sp(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A_dai = [1036.0,1950.0]
	B = 0.044
	r_dai_gt = 0.63
	E_dai = [128000.0,70000.0]
	V_dai = [2.50e3,-0.57e3]

	cond_dry = (A_dai[0] * (1-(B*P)) * np.exp(-((E_dai[0]) + (P*V_dai[0])) / (R_const * T)))

	if (mechanism == None):
		cond = cond_dry
	elif (mechanism == 'polaron'):
		cond = cond_dry
	elif (mechanism == 'dry'):
		cond = cond_dry
	else:
		raise ValueError('Ionic and Proton conduction is not included in electrical conductivity model: Dai2009_DryPyropeGarnet_sp')

	return cond