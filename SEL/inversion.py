#!/usr/bin/env python3

import sys,os
import numpy as np

def _solv_cond_(index, cond_list, object, param, upperlimit, lowerlimit, search_increment, acceptence_threshold, 
	init_guess = None, transition_zone = False, water_solv=False):
	
	
	param_search_array = np.arange(lowerlimit[index], upperlimit[index], search_increment)
	
	if len(param_search_array) == 1:
	
		sol_param = upperlimit[index]
		exec('object.' + param + '[' + str(index) + ']='  + str(upperlimit[index]))
		if water_solv == True:
			if transition_zone == False:
				object.mantle_water_distribute(method = 'index', sol_idx = index)
			else:
				object.transition_zone_water_distribute(method = 'index', sol_idx = index)
		cond_calced = object.calculate_conductivity(method = 'index', sol_idx = index)
		residual = cond_list[index] - cond_calced
	else:

		restart = True
		init_search_increment = np.array(search_increment)
		
		while restart:
			
			restart = False
			
			for j in range(0,len(param_search_array)):
				
				exec('object.' + param + '[' + str(index) + ']='  + str(param_search_array[j]))
				
				if water_solv == True:
					if transition_zone == False:
						object.mantle_water_distribute(method = 'index', sol_idx = index)
					else:
						object.transition_zone_water_distribute(method = 'index', sol_idx = index)
					
				cond_calced = object.calculate_conductivity(method = 'index',sol_idx = index)
				
				residual = cond_list[index] - cond_calced
				
				if abs(residual) < (acceptence_threshold * 1e-2 * cond_list[index]):
					restart = False
					sol_param = param_search_array[j]
					break
					
				else:
					
					if residual < 0.0:
						if len(param_search_array) > 4:
							
							if search_increment <= (init_search_increment * 1e-2 * acceptence_threshold):
								sol_param = lowerlimit[index]
								restart = False
								break
							else:
								search_increment = search_increment / 2.0
								param_search_array = np.arange(param_search_array[j-3], upperlimit[index], search_increment)
								restart = True
								break
						else:					
							if search_increment <= (init_search_increment * 1e-2 * acceptence_threshold):
								sol_param = lowerlimit[index]
								restart = False
								break
							else:
								search_increment = search_increment / 2.0
								param_search_array = np.arange(lowerlimit[index], upperlimit[index], search_increment)
								restart = True
								break
					else:
						if j == len(param_search_array)-1:
							sol_param = param_search_array[-1] #equivalent to upper limit
							restart = False
							break
						else:
							pass
	
	return sol_param, residual
	
def conductivity_solver_single_param(object, cond_list, param_name,
	upper_limit_list, lower_limit_list, search_start, acceptence_threshold, cond_err = None, transition_zone = False, num_cpu = 1):

	index_list = np.array(list(range(0,len(object.T)))) #creating the index array tied to the T array.
	#index_list = 0
	
	if ('water' in param_name) == True:
	
		if param_name == 'bulk_water':
			water_solv = True
			#setting the object.bulk_water as same length as T if that has not done already...
			if len(getattr(object,param_name)) != len(object.T):
				object.set_bulk_water(0.0)
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')
		
		
	elif ('melt' in param_name) == True:
	
		water_solv = True
	
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
							
			process_item_partial = partial(_solv_cond_, cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
			lowerlimit=lower_limit_list , search_increment= search_start, acceptence_threshold = acceptence_threshold, init_guess = 0, transition_zone = transition_zone, water_solv=water_solv)
			
			c = pool.map(process_item_partial, index_list)
			
		c_list = [x[0] for x in c]
		residual_list= [x[1] for x in c]
		
	else:
		
		c_list = np.zeros(len(index_list))
		residual_list = np.zeros(len(index_list))
		
		for idx in range(0,len(index_list)):
			
			c = _solv_cond_(index = index_list[idx], cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
				lowerlimit=lower_limit_list , search_increment= search_start, acceptence_threshold = acceptence_threshold, init_guess = 0, transition_zone = transition_zone, water_solv=water_solv)
			
			c_list[idx] = c[0]
			residual_list[idx] = c[1]
			
	return c_list, residual_list