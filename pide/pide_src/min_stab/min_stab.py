#!/usr/bin/env python3

import numpy as np
from utils.utils import check_type

def Boyd1960_quartz_coesite_trans(T, P):

	type_T = check_type(T)
	
	if type_T == "scalar":
		if T < 973.15:
			idx_coesite = None
		else:
			p_coe = 0.1 * 19.5 * (0.112*(T - 273.15))
			if P < p_coe:
				idx_coesite = None
			else:
				idx_coesite = [0]
	else:
		idx_coesite = []
		for i in range(0,len(T)):
			if T[i] > 973.15:
				p_coe = 0.1 * (19.5 + (0.112*(T[i] - 273.15))) #in GPa
				if P[i] > p_coe:
					idx_coesite.append(i)
					
	return idx_coesite
	
def alpha_beta_quartz(T):
	
	type_T = check_type(T)
	
	if type_T == "scalar":
		
		if T <= 846.5:
			idx_alpha = [0]
		else:
			idx_alpha = None
			
	else:
		idx_alpha = np.where(T<846.5)
		
	return idx_alpha[0]