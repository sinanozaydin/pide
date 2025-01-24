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

def _solv_cond_(index, cond_list, object, param, upperlimit, lowerlimit, search_increment, acceptence_threshold, results_list = None,
	init_guess = None, transition_zone = False, water_solv=False, comp_solv=False, comp_type = None, comp_index = None, low_value_threshold = None,):
	
	"""Simple line-search solver to solve conductivities. This function should not be called directly.	
	"""
	
	if results_list is not None:
		if len(results_list) == 0:
			init_guess = None
		else:
			init_guess = results_list[-1]
			
	param_search_array = np.arange(lowerlimit[index], upperlimit[index] , search_increment)
		
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
		init_restart = True
		
		while restart:
			
			restart = False
			
			if init_guess == None:
				idx_start_search = 0
			else:
				if init_restart == True:
					idx_start_search = np.argmin(np.abs(param_search_array-init_guess))
					init_restart = False
				else:
					idx_start_search = 0
			
			for j in range(idx_start_search,len(param_search_array)):
				
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
						else:
							sol_param = param_search_array[j]

					break
					
				else:
					
					if residual < 0.0:
						if (len(param_search_array) > 4) and (j>=3):
							
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

		if results_list is not None:
			results_list.append(sol_param)
	
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
	
	object.revalue_arrays()
	
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
		
		manager = multiprocessing.Manager()
		shared_results = manager.list()
	
		with multiprocessing.Pool(processes=num_cpu) as pool:
							
			process_item_partial = partial(_solv_cond_, cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
			lowerlimit=lower_limit_list , search_increment= search_start, acceptence_threshold = acceptence_threshold, results_list = shared_results, init_guess = None,
			transition_zone = transition_zone, water_solv=water_solv,comp_solv = comp_solv, comp_type = comp_type, comp_index = comp_index,
			low_value_threshold = low_value_threshold)
			
			c = pool.map(process_item_partial, index_list)
			
		c_list = [x[0] for x in c]
		residual_list= [x[1] for x in c]
							
				
	else:
		
		c_list = np.zeros(len(index_list))
		residual_list = np.zeros(len(index_list))
		
		for idx in range(0,len(index_list)):
			
			if idx > 0:
				init_guess_ = c_list[idx-1]
			else:
				init_guess_ = None
			c = _solv_cond_(index = index_list[idx], cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
				lowerlimit=lower_limit_list , search_increment= search_start, acceptence_threshold = acceptence_threshold, results_list= None, init_guess = init_guess_, transition_zone = transition_zone, 
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
	
#misfit function
def _misfit(cond, cond_external):
	#misfit at log scale
	misf = abs(np.log10(cond) - np.log10(cond_external))
	return misf

# The likelihood function
def _likelihood(cond, cond_external, sigma):
	misf = _misfit(cond, cond_external)
	like = np.exp(-misf**2 / (2 * sigma**2))
	return like, misf

def _solv_MCMC_two_param(index, cond_list, object, initial_params, param_name_1, param_name_2, upper_limits,
	lower_limits, sigma_cond,proposal_stds,n_iter,burning,water_solv, comp_solv, comp_type = ['mineral','mineral'], comp_index = [0,0],
	transition_zone = False, num_cpu = 1):

	#Using Metropolis-Hastings algorithm
	
	frac_1 = False
	frac_2 = False
	
	if 'frac' in param_name_1:
		frac_1 = True
	if 'frac' in param_name_2:
		frac_2 = True
	
	param_1_init, param_2_init = initial_params[index]
	(param_1_max,param_2_max) = upper_limits
	(param_1_min,param_2_min) = lower_limits
	
	current_params = np.array([param_1_init, param_2_init])
	
	#Initial setting of the parameters.
	exec('object.' + param_name_1 + '[' + str(index) + ']='  + str(param_1_init))
	exec('object.' + param_name_2 + '[' + str(index) + ']='  + str(param_2_init))
	
	#Checking if initial water content is entered as other than 0.
	if object.bulk_water[index] > 0.0:
		water_solv = True
	
	if water_solv == True:
	
		if transition_zone == False:
			object.mantle_water_distribute(method = 'index', sol_idx = index)
		else:
			object.transition_zone_water_distribute(method = 'index', sol_idx = index)
	
	#Calculating the initial conductivity
	cond_init = object.calculate_conductivity(method = 'index', sol_idx = index)
	current_likelihood, current_misf = _likelihood(cond_init, cond_list[index], sigma_cond[index])
		
	#empty arrays to fill it up with samples
	samples = []
	misfits = []
	accepted = 0
	
	#loop for monte-carlo
	for _ in range(n_iter):
		#proposing the new parameters
		proposal = current_params + np.random.normal(0, proposal_stds, size=2)
		
		#clipping the proposal distribution, respecting the bounds set up
		proposal[0] = np.clip(proposal[0], param_1_min[index], param_1_max[index])
		proposal[1] = np.clip(proposal[1], param_2_min[index], param_2_max[index])
		
		#setting up the random parameter
		if comp_solv == False:
			exec(f'object.{param_name_1}[{str(index)}] = proposal[0]')
			exec(f'object.{param_name_2}[{str(index)}] = proposal[1]')
		else:
			#Determine which parameter has the
			pass
			
		if water_solv == True:
			if transition_zone == False:
				object.mantle_water_distribute(method = 'index', sol_idx = index)
			else:
				object.transition_zone_water_distribute(method = 'index', sol_idx = index)
				
		proposed_cond = object.calculate_conductivity(method = 'index',sol_idx = index)
		proposed_likelihood, misf = _likelihood(proposed_cond, cond_list[index], sigma_cond[index])
		
		# Calculate acceptance probability
		acceptance_ratio = proposed_likelihood / current_likelihood
		
		if np.random.rand() < acceptance_ratio:
			
			current_params = proposal
			current_likelihood = proposed_likelihood
			misfits.append(misf)
			
			if _ > burning:
				samples.append(current_params)
				accepted += 1
	
	acceptance_rate = accepted / (n_iter - burning)
	
	return np.array(samples), acceptance_rate, misfits

def conductivity_metropolis_hastings_two_param(object, cond_list, initial_params, param_name_1, param_name_2, upper_limits,
	lower_limits, sigma_cond,proposal_stds,n_iter, burning = 0, transition_zone = False,num_cpu = 1):

	"""
	A function to perform stochastic inversion for electrical conductivity with the given two parameters and predefined sets.
	This function utilizes Metropolis-Hastings Markov-Chain Monte Carlo method to produce a posterior distribution of variables.
	
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
	
	if burning >= n_iter:
		
		raise ValueError('Burning samples cannot be larger than the total iteration number (n_iter).')

	if len(cond_list) == len(initial_params) == len(upper_limits[0]) == len(lower_limits[0]) == len(sigma_cond):
		pass
	else:
		raise IndexError('The length of the arrays for each conductivity solution (cond_list) are not same. cond_list, initial_params, upper_limits, lower_limits and sigma_conds has to be the same length.')
	
	min_list = ['quartz_frac', 'plag_frac', 'amp_frac', 'kfelds_frac', 'opx_frac', 'cpx_frac',
		'mica_frac', 'garnet_frac', 'sulphide_frac', 'graphite_frac', 'ol_frac', 'sp_frac', 'rwd_wds_frac',
		'perov_frac', 'mixture_frac', 'other_frac']

	rock_list = ['granite_frac','granulite_frac','sandstone_frac','gneiss_frac','amphibolite_frac',
			'basalt_frac','mud_frac','gabbro_frac','other_rock_frac']
	
	index_list = np.array(list(range(0,len(object.T)))) #creating the index array tied to the T array.
	param_names = [param_name_1, param_name_2]
	
	if any('water' in xx for xx in param_names) == True:
	
		if 'bulk_water' in param_names:
			
			water_solv = True
			comp_solv = False
			comp_type = None
			comp_index = None
			#setting the object.bulk_water as same length as T if that has not done already...
			if len(getattr(object,param_names[param_names.index('bulk_water')])) != len(object.T):
				object.set_bulk_water(0.0)
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')
			
	if ('melt' in xx for xx in param_names) == True:
	
		water_solv = True
		comp_solv = False
		comp_type = None
		comp_index = None
		for ii in range(2):
			if len(getattr(object,param_names[ii])) != len(object.T):
				object.set_parameter(param_names[ii], 0.0)
	
	else:
		water_solv = False
		#setting the object as same length as T if that has not done already...
		
		comp_index = []
		comp_type_list = []
		
		if ('frac' in xx for xx in param_names) == True:
			
			comp_solv = True
			water_solv = True
			
			for ii in range(len(param_names)):
				if param_names[ii] in min_list:
					comp_type = 'mineral'
					comp_type_list.append(comp_type)
					comp_index.append(min_list.index(param_names[ii]))
					
				elif param_names[ii] in rock_list:
					comp_type = 'rock'
					comp_type_list.append(comp_type)
					comp_index.append(rock_list.index(param_names[ii]))
				
				if len(getattr(object,param_names[ii])) != len(object.T):
					object.set_parameter(param_names[ii], 0.0)
				
				else:
					raise NameError('The mineral/rock name you entered is not included as a parameter in pide.')
					
		else:
			for ii in range(len(param_names)):
				if len(getattr(object,param_names[ii])) != len(object.T):
					object.set_parameter(param_names[ii], 0.0)
		
		if (('mineral' in comp_type_list) == True) and (('rock' in comp_type_list) == True):
			raise ValueError('The user cannot enter both rock and mineral as the inversion parameter. Choose only one.')
			
	if num_cpu > 1:
		
		import multiprocessing
		import os
		from functools import partial
		
		max_num_cores = os.cpu_count()
		
		if num_cpu > max_num_cores:
			raise ValueError('There are not enough cpus in the machine to run this action with ' + str(num_cpu) + ' cores.')
			
	if num_cpu > 1:
		
		manager = multiprocessing.Manager()
		shared_results = manager.list()
	
		with multiprocessing.Pool(processes=num_cpu) as pool:
							
			process_item_partial = partial(_solv_MCMC_two_param, object = object, cond_list = cond_list, initial_params = initial_params, param_name_1 = param_name_1, param_name_2= param_name_2,
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond, proposal_stds = proposal_stds , n_iter= n_iter, burning = burning,
			water_solv = water_solv, comp_solv = comp_solv, num_cpu = num_cpu)
			
			c = pool.map(process_item_partial, index_list)
			
		sample_distr = [x[0] for x in c]
		acceptance_rate = [x[1] for x in c]
		misfits = [x[2] for x in c]
		
	else:
	
		sample_distr = []
		acceptance_rate = np.zeros(len(index_list))
		misfits = []
		
		for idx in range(0,len(index_list)):
			
			c = _solv_MCMC_two_param(index = index_list[idx], object = object, cond_list = cond_list, initial_params = initial_params, param_name_1 = param_name_1, param_name_2= param_name_2,
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond, proposal_stds = proposal_stds , n_iter= n_iter, burning = burning,
			water_solv = water_solv, comp_solv = comp_solv, num_cpu = 1)
			
			sample_distr.append(c[0])
			acceptance_rate[idx] = c[1]
			misfits.append(c[2])
			
			
	return sample_distr, acceptance_rate, misfits