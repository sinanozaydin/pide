#!/usr/bin/env python3

def _solv_(index_list, object, param, upperlimit, lowerlimit, search_start, water_solv = False):

	
	
def conductivity_solver_single_param(cond_list, object, temp_array, p_array, param_name,
	upper_limit_list, lower_limit_list, search_start, cond_err = None, num_cpu = 1):

	
	if 'water' in param_name is True:
	
		if param_name == 'bulk_water':
			water_solv = True
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')
	
	else:
		water_solv = False
	
	if num_cpu > 1:
		
		import multiprocessing
		import os
		from functools import partial
		
		max_num_cores = os.cpu_count()
		
		if num_cpu > max_num_cores:
			raise ValueError('There are not enough cpus in the machine to run this action with ' + str(num_cpu) + ' cores.')
			
	if num_cpu > 1:
	
		with multiprocessing.Pool(processes=num_cpu) as pool:
							
			process_item_partial = partial()
			
			c = pool.map(process_item_partial, sliced_material_idx)
	
	else:
	
		water = _solv_()
