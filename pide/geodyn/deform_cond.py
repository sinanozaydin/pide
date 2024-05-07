#!/usr/bin/env python3

import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning

warnings.simplefilter("ignore", OptimizeWarning)

scipy_methods = ['lm', 'trf', 'dogbox']

def check_misfit(x_array, y_array, x_vals, y_vals):
	
	misfit_list = np.zeros(len(x_vals))
	for i in range(0,len(x_vals)):
		misfit_list[i] = y_array[(np.abs(x_array - x_vals[i])).argmin()] - y_vals[i]
	return misfit_list

def plastic_strain_2_conductivity(strain, low_cond, high_cond, low_strain, high_strain, function_method, **kwargs):
	
	"""
	Function to relate accumulated plastic strain to electrical conductivity.
	
	Inputs:
	
	strain: The value to convert.
	
	low_cond: Lower conductivity at lower strain (low_strain) to form the transformation model.
	high_cond: Higher conductivity at higher strain (high_strain) to form the transformation model.
	low_strain: Lower strain value at lower conductivity to form the transformation model.
	high_strain: Higher strain value at higher conductivity to form the transformation model.
	
	function method: Function method to use to build the transformation model. The options are:
	
		1. 'linear'
		2. 'exponential'
		
	Optional parameters:
		
	conductivity_decay_factor: Self-explanatory decay factor (0 to 1) parameter used in exponential method.
	strain_decay_factor: Self-explanatory decay factor (0 to 1) parameter used in exponential method.
	
	These decay factors controls how the function decays between lower values and higher values.
	
	return_all_params: Function returns all relevant parameters instead of just conductivity.
	"""
	
	
	if function_method == 'exponential':
		cond_decay_modifier = -1
	else:
		cond_decay_modifier = 1
		
	conductivity_decay_factor = cond_decay_modifier * kwargs.pop('conductivity_decay_factor', 0.0)
	strain_decay_factor = kwargs.pop('strain_decay_factor', 0.0)
	return_all_params = kwargs.pop('return_all_params', False)
	strain_percolation_threshold = kwargs.pop('strain_percolation_threshold', None)
	
	#Error values as kwargs:
	low_cond_err = kwargs.pop('low_cond_err', 10)
	mid_cond_err = kwargs.pop('mid_cond_err', 10)
	high_cond_err = kwargs.pop('high_cond_err', 10)
	
	strain_func_build = np.logspace(np.log10(low_strain),np.log10(high_strain),100)

	def setup_data(low_strain, high_strain, low_cond, high_cond, strain_decay_factor, conductivity_decay_factor, low_cond_err = low_cond_err, 
					high_cond_err = high_cond_err, mid_cond_err = mid_cond_err):
	
		strains = np.array([low_strain, high_strain])
		conds = np.array([low_cond, high_cond])
		
		strains_log = np.log10(strains)
		conds_log = np.log10(conds)
		
		mid_strains_log = ((strains_log[1] - strains_log[0]) / 2.0) + strains_log[0]
		mid_conds_log = ((conds_log[1] - conds_log[0]) / 2.0) + conds_log[0]
	
		if (strain_decay_factor < -1) or (strain_decay_factor >= 1.0):
	
			raise ValueError('The value entered for the x-factor is not valid. It has to be in between 0 to 1')
		
		if (conductivity_decay_factor < -1.0) or (conductivity_decay_factor >= 1.0):
	
			raise ValueError('The value entered for the y-factor is not valid. It has to be in between 0 to 1')
		
		else:
						
			mid_conds_log = ((mid_conds_log - conds_log[0]) * conductivity_decay_factor) + mid_conds_log
			mid_strains_log = ((mid_strains_log - strains_log[0]) * strain_decay_factor) + mid_strains_log
			strains = 10**np.array([strains_log[0],mid_strains_log,strains_log[1]])
			conds = 10**np.array([conds_log[0],mid_conds_log,conds_log[1]])
			cond_err = np.array([conds[0]* (low_cond_err / 1e2), conds[1] * (mid_cond_err / 1e2), conds[2] * (high_cond_err / 1e2)])			
	
		return strains, conds, cond_err
		
	strains, conds, cond_err = setup_data(low_strain, high_strain, low_cond, high_cond, strain_decay_factor, conductivity_decay_factor)
	
	if function_method == 'linear':
	
		f = linregress(strains,conds)
	
		conds_calc = f.intercept + (f.slope*strain)
		
		return conds_calc, misfit
		
	elif function_method == 'exponential':
		
		#setting up run iteration number to check what to do for convergence of the function
		run_idx = 0
	
		#exponential function definition
		def func(x,a,b,c):
		
			return (a * np.exp(b*x)) + c
								
		def fit(method,init_start, run_idx, strains, conds, strain_decay_factor, conductivity_decay_factor, cond_err, func_method = 'exponential'):
			try:
				if init_start == True:
					run_idx = 0
				else:
					run_idx = run_idx
					
				if func_method == 'exponential':
					
					try:
						params, params_cov = curve_fit(func, strains, conds, method = method, sigma = cond_err)
					except OptimizeWarning as e:
						pass
					cond_calced = func(strain_func_build, params[0], params[1], params[2])
					misfit = check_misfit(np.log10(strains),np.log10(cond_calced),np.log10(strains),np.log10(conds)) #log misfit
					
				elif func_method == 'linear':
					
					f = linregress(strains,conds)
					cond_calced = f.intercept + (f.slope*strain_func_build)
					misfit = check_misfit(np.log10(strains),np.log10(cond_calced),np.log10(strains),np.log10(conds)) #log misfit
				
				if np.any(misfit > 0.5):
				
					raise RuntimeError
					
				else:
					
					if func_method == 'exponential':
						cond_calced = func(strain, params[0], params[1], params[2])
					elif func_method == 'linear':
						cond_calced = f.intercept + (f.slope*strain)
					if strain_percolation_threshold != None:
						cond_calced[np.where(strain>strain_percolation_threshold)[0]] = conds[-1]
					
			except RuntimeError:
				run_idx = run_idx + 1
				if run_idx < -5:
					cond_calced, strain_decay_factor, conductivity_decay_factor, misfit = fit(method = scipy_methods[run_idx], strains = strains, conds = conds, 
					init_start=False, run_idx = run_idx, strain_decay_factor=strain_decay_factor, conductivity_decay_factor=conductivity_decay_factor, cond_err = cond_err)
				elif run_idx >= 1:
					#modifying the decay factor slightly to fit making more exponential 
					if strain_decay_factor < 0.75:
						strain_decay_factor = strain_decay_factor * 1.05
					else:
						strain_decay_factor = strain_decay_factor * 0.95
						
					if conductivity_decay_factor < 0.75:
						conductivity_decay_factor = conductivity_decay_factor * 1.05
					else:
						conductivity_decay_factor = conductivity_decay_factor * 0.95
					
					strains, conds, cond_err = setup_data(low_strain, high_strain, low_cond, high_cond, strain_decay_factor, conductivity_decay_factor)
					
					if run_idx < 10:
						cond_calced,strain_decay_factor,conductivity_decay_factor, misfit = fit(method = scipy_methods[0], init_start = False, run_idx = run_idx,
						strains = strains, conds = conds, strain_decay_factor=strain_decay_factor, conductivity_decay_factor=conductivity_decay_factor, cond_err = cond_err)
					else:
						cond_calced,strain_decay_factor,conductivity_decay_factor, misfit = fit(method = scipy_methods[0], init_start = False, run_idx = run_idx, 
						strains = strains, conds = conds, strain_decay_factor=strain_decay_factor, conductivity_decay_factor=conductivity_decay_factor, cond_err = cond_err, func_method= 'linear')
					
				else:
					raise RuntimeError('The optimal parameters for the fit cannot be found. Try adjusting x-factor and y-factor.')
			
			return cond_calced, strain_decay_factor, conductivity_decay_factor, misfit
		
		
		cond_calced,strain_decay_factor,conductivity_decay_factor,misfit = fit(method = scipy_methods[0], init_start = True, run_idx = 0,
		strains = strains, conds = conds, strain_decay_factor= strain_decay_factor, conductivity_decay_factor=conductivity_decay_factor, cond_err = cond_err)
		
		if return_all_params == False:
			return cond_calced
		elif return_all_params == True:
			return cond_calced, strain_decay_factor, conductivity_decay_factor, misfit[1]
			