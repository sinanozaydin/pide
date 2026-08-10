#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

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