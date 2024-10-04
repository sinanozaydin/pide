#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Duba1978_Archie_param1_porosity(T,P,water, param1, fo2 = None, fo2_ref = None, method = None):

	#Taking from porosity fitting to Duba 1978 equation 5
	#param1 is porosity in volumetric percentage...
	
	if np.any(param1 == 0.0):
		raise ValueError('param1 for sandstone has to be set to a value for this function. It denotes to porosity volumetric percentage in sandstone.')
	
	cond = 10.0**((1.214*np.log10(param1)) - 2.72)
	
	return cond