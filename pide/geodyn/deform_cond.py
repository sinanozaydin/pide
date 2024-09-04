#!/usr/bin/env python3

import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning
import matplotlib.pyplot as plt

warnings.simplefilter("ignore", OptimizeWarning)

scipy_methods = ['lm', 'trf', 'dogbox']

def check_misfit(x_array, y_array, x_vals, y_vals, y_err):
	
	misfit_list = np.zeros(len(x_vals))
	for i in range(0,len(x_vals)):
		misfit_list[i] = y_array[(np.abs(x_array - x_vals[i])).argmin()] - y_vals[i]
	misfitsum = np.sum((misfit_list/y_err)**2.0)
	rms = np.sqrt((1.0 / len(x_vals)) * misfitsum)

	return rms, misfit_list

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
		3. 'polynomial'
		
	Optional parameters:
		
	conductivity_decay_factor: Self-explanatory decay factor (0 to 1) parameter used in exponential method.
	strain_decay_factor: Self-explanatory decay factor (0 to 1) parameter used in exponential method.
	
	These decay factors controls how the function decays between lower values and higher values.
	
	return_all_params: Function returns all relevant parameters instead of just conductivity.
	"""
	
			
	conductivity_decay_factor = kwargs.pop('conductivity_decay_factor', 0.7)
	conductivity_decay_factor_2 = kwargs.pop('conductivity_decay_factor_2', 0.2)
	strain_decay_factor = kwargs.pop('strain_decay_factor', 0.3)
	return_all_params = kwargs.pop('return_all_params', False)
	strain_percolation_threshold = kwargs.pop('strain_percolation_threshold', None)
	index_plot = kwargs.pop('index_plot', 1)
	
	#Error values as kwargs:
	cond_err = kwargs.pop('cond_err', 0.1)
	
	strain_func_build = np.log10(np.logspace(np.log10(low_strain),np.log10(high_strain),50))

	def setup_data(low_strain, high_strain, low_cond, high_cond, strain_decay_factor, conductivity_decay_factor, conductivity_decay_factor_2, cond_err = cond_err):
					
		if (strain_decay_factor < -0.9) or (strain_decay_factor >= 0.9):
	
			raise ValueError('The value entered for the strain_decay_factor is not valid. It has to be in between 0 to 0.9')
		
		if (conductivity_decay_factor < -1.0) or (conductivity_decay_factor >= 1.0):
	
			raise ValueError('The value entered for the conductivity_decay_factor is not valid. It has to be in between 0 to 1')
			
		if (conductivity_decay_factor_2 < -1.0) or (conductivity_decay_factor_2 >= 1.0):
	
			raise ValueError('The value entered for the conductivity_decay_factor_2 is not valid. It has to be in between 0 to 1')
			
		if abs(conductivity_decay_factor + conductivity_decay_factor_2) > 1.0:
		
			raise ValueError('The sum value of conductivity_decay_factor(s) cannot be over 1.')
	
		strains = np.array([low_strain, high_strain])
		conds = np.array([low_cond, high_cond])
		
		strains_log = np.log10(strains)
		conds_log = np.log10(conds)
		
		#finding middle point
		mid_strains_log = ((strains_log[1] - strains_log[0]) / 2.0) + strains_log[0]
		mid_conds_log = ((conds_log[1] - conds_log[0]) / 2.0) + conds_log[0]
				
		first_strains_log = mid_strains_log - ((mid_strains_log - strains_log[0]) * strain_decay_factor)
		second_strains_log = mid_strains_log + ((mid_strains_log - strains_log[0]) * strain_decay_factor)
		
		center_conds_log = mid_conds_log - ((mid_conds_log - conds_log[0]) * conductivity_decay_factor)
		first_conds_log = center_conds_log - ((mid_conds_log - conds_log[0]) * conductivity_decay_factor_2)
		second_conds_log = center_conds_log + ((mid_conds_log - conds_log[0]) * conductivity_decay_factor_2)
		
		strains = np.array([strains_log[0],first_strains_log,second_strains_log,strains_log[1]])
		conds = np.array([conds_log[0],first_conds_log,second_conds_log,conds_log[1]])
	
		conds_err = abs(conds * cond_err)
		
		return strains, conds, conds_err
		
	strains, conds, cond_err = setup_data(low_strain, high_strain, low_cond, high_cond, strain_decay_factor,
	conductivity_decay_factor, conductivity_decay_factor_2,cond_err)
	
	if function_method == 'linear':
		
		f = linregress([strains[0],strains[-1]],[conds[0],conds[-1]])
	
		cond_calced = 10.0**(f.intercept + (f.slope*np.log10(strain)))
		rms = 0.0
		
		if return_all_params == False:
			return cond_calced
		elif return_all_params == True:
			return cond_calced, rms
		
	elif function_method in ['exponential', 'polynomial']:
		
		#setting up run iteration number to check what to do for convergence of the function
		run_idx = 0
	
		#exponential function definition
		def expfunc(x,a,b,c):
		
			return (a * np.exp(b*x)) + c
			
		def polyfunc(x,a,b,c,d):
		
			return a + (b*x) + (c*x**2.0) + (d*x**3.0)
								
		def fit(method,init_start, run_idx, strains, conds, cond_err, func_method = 'exponential'):
			
			#Subfunction that fits the exponantial curve recursively.
		
			if init_start == True:
				run_idx = 0
			else:
				run_idx = run_idx
				
			if func_method in ['exponential', 'polynomial']:
				
				try:
					if func_method == 'exponential':
						params, params_cov = curve_fit(expfunc, strains, conds, method = method, sigma = cond_err)
					elif func_method == 'polynomial':
						params, params_cov = curve_fit(polyfunc, strains, conds, method = method, sigma = cond_err)
				except OptimizeWarning as e:
					pass
				if func_method == 'exponential':
					cond_calced = expfunc(strain_func_build, params[0], params[1], params[2])
					rms, misfit_list = check_misfit(strain_func_build,cond_calced,strains,conds,cond_err)
				elif func_method == 'polynomial':
					cond_calced = polyfunc(strain_func_build, params[0], params[1], params[2], params[3])
					rms, misfit_list = check_misfit(strain_func_build,cond_calced,strains,conds,cond_err)
					
			else: 
				
				if func_method == 'linear':
					f = linregress(strains,conds)
					cond_calced = f.intercept + (f.slope*strain_func_build)
					rms, misfit_list = check_misfit(np.log10(strains),np.log10(cond_calced),np.log10(strains),np.log10(conds),np.log10(cond_err)) #log misfit
					
				else:
					
					raise ValueError('The function method can only be one of those string values: exponential, logarithmic, linear.')
			
			"""
			if rms > 1.5:
				print('RMS values higher than 1.5 is detected, this does not mean an unreasonable fit is achived. Mapping RMS values may help understand where this problem	is occuring.')			
			"""
			
			if func_method == 'exponential':
				cond_calced = 10.0**expfunc(np.log10(strain), params[0], params[1], params[2])
				
			elif func_method == 'polynomial':
				cond_calced = 10.0**polyfunc(np.log10(strain), params[0], params[1], params[2], params[3])
							
			return cond_calced, rms
		
		cond_calced, rms = fit(method = 'lm', init_start = True, run_idx = 0,
		strains = strains, conds = conds, cond_err = cond_err, func_method=function_method)
		
		if return_all_params == False:
			return cond_calced
		elif return_all_params == True:
			return cond_calced, rms