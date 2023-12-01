#!/usr/bin/env python3

import sys,os
import numpy as np

def _solv_(index, cond_list, object, param, upperlimit, lowerlimit, search_start, init_guess, water_solv=False):
					
	if water_solv == True:
	
		exec('object.' + param + '[' + str(index) + ']='  + str(lowerlimit[index]))
		object.mantle_water_distribute(method = 'index')
		
	else:
		object.set_parameter(lowerlimit[index])
		
	cond_calced = object.calculate_conductivity(method = 'index')
	
	residual = cond_list[index] - cond_calced[index]
	
	param_search_array = np.arange(lowerlimit[index], upperlimit[index], search_start)
	
	
	
def conductivity_solver_single_param(object, cond_list, param_name,
	upper_limit_list, lower_limit_list, search_start, cond_err = None, num_cpu = 1):

	# index_list = np.array(list(range(0,len(object.T)))) #creating the index array tied to the T array.
	index_list = 0
	
	if ('water' in param_name) == True:
	
		if param_name == 'bulk_water':
			water_solv = True
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')
		
		#setting the object.bulk_water as same length as T if that has not done already...
		if len(getattr(object,param_name)) != len(object.T):
			object.set_bulk_water(0.0)
	
	else:
		water_solv = False
		#setting the object as same length as T if that has not done already...
		if len(getattr(object,param_name)) != len(object.T):
			object.set_parameter(param_name, 0.0)
	
	if num_cpu > 1:
		
		import multiprocessing
		import os
		from functools import partial
		
		max_num_cores = os.cpu_count()
		
		if num_cpu > max_num_cores:
			raise ValueError('There are not enough cpus in the machine to run this action with ' + str(num_cpu) + ' cores.')
			
	if num_cpu > 1:
	
		with multiprocessing.Pool(processes=num_cpu) as pool:
							
			process_item_partial = partial(_solv_, cond_list, temp_array, pres_array, object, param_name, upper_limit_list, lower_limit_list, search_start)
			
			c = pool.map(process_item_partial, index_list)
	
	else:
	
		water = _solv_(index = index_list, cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
		lowerlimit=lower_limit_list , search_start= search_start, init_guess = 0, water_solv=water_solv)
