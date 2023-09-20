#!/usr/bin/env python3

import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit

scipy_methods = ['lm', 'trf', 'dogbox']

def check_misfit(x_array, y_array, x_vals, y_vals):
	
	misfit_list = np.zeros(len(x_vals))
	for i in range(0,len(x_vals)):
		misfit_list[i] = y_array[np.where(x_array == x_vals[i])[0][0]] - y_vals[i]
	return misfit_list

def plastic_strain_2_conductivity(strain, low_cond, high_cond, low_strain, high_strain, function_method, **kwargs):

	conductivity_decay_factor = -kwargs.pop('conductivity_decay_factor', 0.0)
	strain_decay_factor = kwargs.pop('strain_decay_factor', 0.0)

	strains = np.array([low_strain, high_strain])
	conds = np.array([low_cond, high_cond])
	
	strains_log = np.log10(strains)
	conds_log = np.log10(conds)
	strain_log = np.log10(strain)
	
	if function_method == 'linear':
	
		f = linregress(strains,conds)
	
		conds_calc = f.intercept + (f.slope*strain)
		
		return conds_calc, misfit
		
	elif function_method == 'exponential':
		
		#
		run_idx = 0
		
		mid_strains_log = ((strains_log[1] - strains_log[0]) / 2.0) + strains_log[0]
		mid_conds_log = ((conds_log[1] - conds_log[0]) / 2.0) + conds_log[0]
	
		if (strain_decay_factor < 0) or (strain_decay_factor >= 1.0):
	
			raise ValueError('The value entered for the x-factor is not valid. It has to be in between 0 to 1')
		
		if (conductivity_decay_factor < -1.0) or (conductivity_decay_factor >= 1.0):
	
			raise ValueError('The value entered for the y-factor is not valid. It has to be in between 0 to 1')
		
		else:
	
			mid_conds_log = ((mid_conds_log - conds_log[0]) * conductivity_decay_factor) + mid_conds_log
			mid_strains_log = ((mid_strains_log - strains_log[0]) * strain_decay_factor) + mid_strains_log
			strains = 10**np.array([strains_log[0],mid_strains_log,strains_log[1]])
			conds = 10**np.array([conds_log[0],mid_conds_log,conds_log[1]])
		
		#exponential function definition
		def func(x,a,b,c):
			return (a * np.exp(b*x)) + c
		
		def linear_func(x,c,d):
		
			return (c*x) + d
						
		def fit(method,init_start, run_idx, strain_decay_factor, conductivity_decay_factor):
			try:
				if init_start == True:
					run_idx = 0
				else:
					run_idx = run_idx
	
				params, params_cov = curve_fit(func, strains, conds, method = method)
				cond_calced = func(strain, params[0], params[1], params[2])
				misfit = check_misfit(np.log10(strains),np.log10(cond_calced),np.log10(strains),np.log10(conds)) #log misfit
		
				if np.any(misfit > 0.5):
					
					raise RuntimeError
					
			except RuntimeError:
				run_idx = run_idx + 1
				if run_idx < 3:
					fit(method = scipy_methods[run_idx], init_start=False, run_idx = run_idx, strain_decay_factor=strain_decay_factor, conductivity_decay_factor=conductivity_decay_factor)
				elif run_idx >= 3:
					#changing 
					strain_decay_factor = strain_decay_factor * 0.9
					conductivity_decay_factor = conductivity_decay_factor * 0.9
					run_idx = 0
					fit(method = scipy_methods[run_idx], init_start = False, run_idx = run_idx, strain_decay_factor=strain_decay_factor, conductivity_decay_factor=conductivity_decay_factor)
				else:
					raise RuntimeError('The optimal parameters for the fit cannot be found. Try adjusting x-factor and y-factor.')
	
			
			
			return cond_calced
			
		cond_calced = fit(method = scipy_methods[0], init_start = True, run_idx = 0, strain_decay_factor= strain_decay_factor, conductivity_decay_factor=conductivity_decay_factor)
		
		return cond_calced