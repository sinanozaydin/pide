#!/usr/bin/env python3

import numpy as np
R_const = 8.3144621

def Liu2020_OpxSol(T,P,depth,h2o_fug,o2_fug,fe_opx,al_opx,method):
	
	h2o_fug = h2o_fug * 1e5 #converting to bars
	
	n = 0.5
	dv = 6.6
	E = 22.5
			
	max_opx_water = (h2o_fug**n) * np.exp(-(E + P*dv) / (R_const * T))
			
	return max_opx_water
	
def Liu2020_OpxSol_n1(T,P,depth,h2o_fug,o2_fug,fe_opx,al_opx,method):
	
	h2o_fug = h2o_fug * 1e5 #converting to bars
	
	n = 1.0
	dv = 16.1e3
	E = 10100.0
			
	max_opx_water = (h2o_fug**n) * np.exp(-(E + P*dv) / (R_const * T))
			
	return max_opx_water
