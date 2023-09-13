#!/usr/bin/env python3

import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit

def plastic_strain_2_conductivity(strain, low_cond, high_cond, low_strain, high_strain, function_method, **kwargs):

	strains = np.array([low_strain, high_strain])
	conds = np.array([low_cond, high_cond])

	if function_method == 'linear':
	
		f = linregress(strains,conds)
	
		conds_calc = f.intercept + (f.slope*strain)
		
		return conds_calc
		
	elif function_method == 'exponential':
	
		conds_calc = 0
		
		return conds_calc