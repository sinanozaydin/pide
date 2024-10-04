#!/usr/bin/env python3

import numpy as np

def _comp_adjust_(_comp_list, comp_alien, comp_old,final = False):

	"""A method to adjust composition of one mineral/rock without considering the replacement weights
	"""
	
	if final == False:
		ratio = (comp_alien - comp_old) / (np.sum(_comp_list) - comp_old)
	else:
		ratio = (comp_alien - comp_old) / (np.sum(_comp_list,axis = 0) - comp_old)
	comp_list = _comp_list - (_comp_list * ratio)
	
	return comp_list

def _solv_cond_(index, cond_list, object, param, upperlimit, lowerlimit, search_increment, acceptence_threshold, 
	init_guess = None, transition_zone = False, water_solv=False, comp_solv=False, comp_type = None, comp_index = None, low_value_threshold = None):
	
	"""Simple grid-search solver to solve conductivities. This function should not be called directly.	
	"""
		
	param_search_array = np.arange(lowerlimit[index], upperlimit[index], search_increment)
	
	if len(param_search_array) == 1:
	
		sol_param = upperlimit[index]
		if comp_solv == True:
			if comp_type == 'mineral':
				_comp_list = [object.quartz_frac[index], object.plag_frac[index], object.amp_frac[index], object.kfelds_frac[index], object.opx_frac[index], object.cpx_frac[index],
				object.mica_frac[index], object.garnet_frac[index], object.sulphide_frac[index], object.graphite_frac[index], object.ol_frac[index], object.sp_frac[index], object.rwd_wds_frac[index],
				object.perov_frac[index], object.mixture_frac[index], object.other_frac[index]]
				comp_old = _comp_list[comp_index]
			elif comp_type == 'rock':
				_comp_list = [object.granite_frac[index],object.granulite_frac[index],object.sandstone_frac[index],object.gneiss_frac[index],object.amphibolite_frac[index],
				object.basalt_frac[index],object.mud_frac[index],object.gabbro_frac[index],object.other_rock_frac[index]]
							
			comp_list = _comp_adjust_(np.array(_comp_list), param_search_array[j], comp_old)
			
			for idx_t in range(len(_comp_list)):
				
				if comp_type == 'mineral':
					object.mineral_frac_list[idx_t][index] = comp_list[idx_t]
				elif comp_type == 'rock':
					object.rock_frac_list[idx_t][index] = comp_list[idx_t]
	
			exec('object.' + param + '[' + str(index) + ']='  + str(upperlimit[index]))
			
			if object.bulk_water[index] > 0.0:
				water_solv = True
					
		else:
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
			
				if comp_solv == True:
				
					if comp_type == 'mineral':
						_comp_list = [object.quartz_frac[index], object.plag_frac[index], object.amp_frac[index], object.kfelds_frac[index], object.opx_frac[index], object.cpx_frac[index],
						object.mica_frac[index], object.garnet_frac[index], object.sulphide_frac[index], object.graphite_frac[index], object.ol_frac[index], object.sp_frac[index], object.rwd_wds_frac[index],
						object.perov_frac[index], object.mixture_frac[index], object.other_frac[index]]
						comp_old = _comp_list[comp_index]
					elif comp_type == 'rock':
						_comp_list = [object.granite_frac[index],object.granulite_frac[index],object.sandstone_frac[index],object.gneiss_frac[index],object.amphibolite_frac[index],
						object.basalt_frac[index],object.mud_frac[index],object.gabbro_frac[index],object.other_rock_frac[index]]
										
					comp_list = _comp_adjust_(np.array(_comp_list), param_search_array[j], comp_old)
					
					for idx_t in range(len(_comp_list)):
						
						if comp_type == 'mineral':
							object.mineral_frac_list[idx_t][index] = comp_list[idx_t]
						elif comp_type == 'mineral':
							object.rock_frac_list[idx_t][index] = comp_list[idx_t]
						
					exec('object.' + param + '[' + str(index) + ']='  + str(param_search_array[j]))
					
					if object.bulk_water[index] > 0.0:
						water_solv = True
					else:
						water_solv = False
						
				else:
				
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
					if low_value_threshold is None:
						sol_param = param_search_array[j]
					else:
						if param_search_array[j] < low_value_threshold:
						
							sol_param = 0.0
							
							if comp_solv == True:
							
								comp_list = _comp_adjust_(np.array(_comp_list), 0.0 , comp_old)
								
								for idx_t in range(len(_comp_list)):
						
									if comp_type == 'mineral':
										object.mineral_frac_list[idx_t][index] = comp_list[idx_t]
									elif comp_type == 'mineral':
										object.rock_frac_list[idx_t][index] = comp_list[idx_t]
									
								exec('object.' + param + '[' + str(index) + ']='  + str(param_search_array[j]))
						else:
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
	upper_limit_list, lower_limit_list, search_start, acceptence_threshold, cond_err = None, transition_zone = False, num_cpu = 1,**kwargs):

	"""
	A function to fit conductivity value with a single parameter with simple line-search algorithm.
	
	Input:
	object: object - pide object that is going to be used for the inversion calculations.
	array: cond_list - conductivity array list used for inversion || in S/m
	str: param_name - parameter name to invert for. This can be any parameter that is included in pide.object.
	array: upper_limit_list - upper limit value for the search space for the given parameter.
	array: lower_limit_list - lower limit value for the search space for the given parameter.
	float: search_start - initial search length used in the line search.
	float: acceptance_threshold - acceptance value to stop inversion process.
	array: cond_err - error floors to add to the inversion.
	bool: transition_zone - boolean value to indicate transition zone water distribution functions are going to be used.
	int: num_cpu - number of cpu to compute the inversion.
	float: low_value_threshold - threshold value of parameter to revert to zero for the solution. 
	"""
	
	min_list = ['quartz_frac', 'plag_frac', 'amp_frac', 'kfelds_frac', 'opx_frac', 'cpx_frac',
		'mica_frac', 'garnet_frac', 'sulphide_frac', 'graphite_frac', 'ol_frac', 'sp_frac', 'rwd_wds_frac',
		'perov_frac', 'mixture_frac', 'other_frac']

	rock_list = ['granite_frac','granulite_frac','sandstone_frac','gneiss_frac','amphibolite_frac',
			'basalt_frac','mud_frac','gabbro_frac','other_rock_frac']

	index_list = np.array(list(range(0,len(object.T)))) #creating the index array tied to the T array.
	
	low_value_threshold = kwargs.pop('low_value_threshold', None)
	
	if ('water' in param_name) == True:
	
		if param_name == 'bulk_water':
			water_solv = True
			comp_solv = False
			comp_type = None
			comp_index = None
			#setting the object.bulk_water as same length as T if that has not done already...
			if len(getattr(object,param_name)) != len(object.T):
				object.set_bulk_water(0.0)
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')
			
	elif ('melt' in param_name) == True:
	
		water_solv = True
		comp_solv = False
		comp_type = None
		comp_index = None
		if len(getattr(object,param_name)) != len(object.T):
			object.set_parameter(0.0)
	
	else:
		water_solv = False
		#setting the object as same length as T if that has not done already...
		if len(getattr(object,param_name)) != len(object.T):
			object.set_parameter(param_name, 0.0)
	
		if ('frac' in param_name) == True:
			
			comp_solv = True
			
			if param_name in min_list:
				comp_type = 'mineral'
				comp_index = min_list.index(param_name)
			elif param_name in rock_list:
				comp_type = 'rock'
				comp_index = rock_list.index(param_name)
			else:
				raise NameError('The mineral/rock name you entered is not included as a parameter in pide.')
	
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
			lowerlimit=lower_limit_list , search_increment= search_start, acceptence_threshold = acceptence_threshold, init_guess = 0, transition_zone = transition_zone,
			water_solv=water_solv,comp_solv = comp_solv, comp_type = comp_type, comp_index = comp_index, low_value_threshold = low_value_threshold)
			
			c = pool.map(process_item_partial, index_list)
			
		c_list = [x[0] for x in c]
		residual_list= [x[1] for x in c]
							
				
	else:
		
		c_list = np.zeros(len(index_list))
		residual_list = np.zeros(len(index_list))
		
		for idx in range(0,len(index_list)):
			
			c = _solv_cond_(index = index_list[idx], cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
				lowerlimit=lower_limit_list , search_increment= search_start, acceptence_threshold = acceptence_threshold, init_guess = 0, transition_zone = transition_zone, 
				water_solv=water_solv, comp_solv = comp_solv, comp_type = comp_type, comp_index = comp_index, low_value_threshold = low_value_threshold)
			
			c_list[idx] = c[0]
			residual_list[idx] = c[1]
			
	#Block to assign the values
	if comp_solv == True:
		
		if param_name in min_list:
			exec(f'comp_list = _comp_adjust_(object.mineral_frac_list, c_list, object.{param_name}, final = True)')
			
			for idx_set in range(len(object.mineral_frac_list)):
			
				if idx_set == comp_index:
					eval(f"object.set_parameter('{param_name}',c_list)")
				else:
					eval(f"object.set_parameter('{min_list[idx_set]}',comp_list[idx_set])")
					
		elif param_name in rock_list:
			
			exec(f'comp_list = _comp_adjust_(object.rock_frac_list, c_list, object.{param_name}, final = True)')
			
			for idx_set in range(len(object.rock_frac_list)):
			
				if idx_set == comp_index:
					eval(f"object.set_parameter('{param_name}',c_list)")
				else:
					eval(f"object.set_parameter('{rock_list[idx_set]}',comp_list[idx_set])")
			
	return c_list, residual_list
	