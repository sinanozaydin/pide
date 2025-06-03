#!/usr/bin/env python3

import numpy as np
from .utils.utils import text_color, check_type

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
	init_guess = None, transition_zone = False, water_solv=False, comp_solv=False, melt_solv=False, comp_type = None, comp_index = None, low_value_threshold = None,
	sfd = False):

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

			print(text_color.RED + 'WARNING: There is no search array can be created with the given lower/upper limit and search increment for composition. Try to change these parameters to have a more reliable solution.' + text_color.END)

			if comp_type == 'mineral':
				_comp_list = [object.quartz_frac[index], object.plag_frac[index], object.amp_frac[index], object.kfelds_frac[index], object.opx_frac[index], object.cpx_frac[index],
				object.mica_frac[index], object.garnet_frac[index], object.sulphide_frac[index], object.graphite_frac[index], object.ol_frac[index], object.sp_frac[index], object.rwd_wds_frac[index],
				object.perov_frac[index], object.mixture_frac[index], object.other_frac[index]]
				comp_old = _comp_list[comp_index]
			elif comp_type == 'rock':
				_comp_list = [object.granite_frac[index],object.granulite_frac[index],object.sandstone_frac[index],object.gneiss_frac[index],object.amphibolite_frac[index],
				object.basalt_frac[index],object.mud_frac[index],object.gabbro_frac[index],object.other_rock_frac[index]]

			comp_list = _comp_adjust_(np.array(_comp_list), upperlimit[index], comp_old)

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
		cond_calced = object.calculate_conductivity(method = 'index', sol_idx = index, sfd = sfd)
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

				cond_calced = object.calculate_conductivity(method = 'index',sol_idx = index, sfd = sfd)
				
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
								#sol_param = lowerlimit[index]
								sol_param = param_search_array[0]
								restart = False
								break
							else:
								search_increment = search_increment / 2.0
								param_search_array = np.arange(param_search_array[j-3], upperlimit[index], search_increment)
								restart = True
								break
						else:
							if search_increment <= (init_search_increment * 1e-2 * acceptence_threshold):
								sol_param = param_search_array[0]
								#sol_param = lowerlimit[index]
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
	upper_limit_list, lower_limit_list, search_start, acceptence_threshold, cond_err = None, transition_zone = False, simplify_fluid_density = False,
	num_cpu = 1,**kwargs):

	"""
	Fit a single parameter to conductivity data using a simple line-search algorithm.

	Parameters
	----------
	object : object
		An instance of the pide object used for inversion calculations.
	cond_list : array-like
		List or array of observed conductivity values [S/m].
	param_name : str
		Name of the parameter to invert. Must be an attribute of the input object.
	upper_limit_list : array-like
		Upper bounds for the search space of the parameter.
	lower_limit_list : array-like
		Lower bounds for the search space of the parameter.
	search_start : float
		Initial step size for the line search.
	acceptence_threshold : float
		Convergence threshold in % of the value entered in cond_list; the search stops when improvement is below this value.
	cond_err : array-like, optional
		Error floors to apply to the conductivity values [S/m] (default is None).
	transition_zone : bool, optional
		Whether to use transition zone water distribution functions (default is False).
	simplify_fluid_density : bool, optional
		If True, simplify fluid density calculation (default is False).
	num_cpu : int, optional
		Number of CPU cores to use in the inversion (default is 1).
	low_value_threshold : float, optional
		Threshold below which parameter values are treated as zero (default is None).
	melt_solv : bool, optional
		If True, include melt solubility in the model (default is False).
		
	Returns:
	-------
	
	c_list: array-like
		solution to the inversion for the chosen parameter
	residuals: array-like
		residuals from the solution.
		
	Examples:
	
	Example for solving bulk water content:
	
	conductivity_solver_single_param(object=object, cond_list = [0.1,0.1], param_name = 'bulk_water',
	upper_limit_list = [1000,1000], lower_limit_list = [0,0],
	search_start = 30, acceptence_threshold = 1, cond_err = None, transition_zone = False, simplify_fluid_density = False,
	num_cpu = 5, melt_solv = 0, low_value_threshold = 10)
	"""

	min_list = ['quartz_frac', 'plag_frac', 'amp_frac', 'kfelds_frac', 'opx_frac', 'cpx_frac',
		'mica_frac', 'garnet_frac', 'sulphide_frac', 'graphite_frac', 'ol_frac', 'sp_frac', 'rwd_wds_frac',
		'perov_frac', 'mixture_frac', 'other_frac']

	rock_list = ['granite_frac','granulite_frac','sandstone_frac','gneiss_frac','amphibolite_frac',
			'basalt_frac','mud_frac','gabbro_frac','other_rock_frac']

	index_list = np.array(list(range(0,len(object.T)))) #creating the index array tied to the T array.

	low_value_threshold = kwargs.pop('low_value_threshold', None)
	melt_solv = kwargs.pop('melt_solv', False)

	object.revalue_arrays()
	
	if check_type(cond_list) != 'array':
		raise KeyError('The value entered for cond_list has to be a list or a numpy array matching the length of temperature array.')

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
		melt_solv = True
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
				
		else:
			
			comp_solv = False
			comp_type = None
			comp_index = None

	if num_cpu > 1:

		import multiprocessing
		import os
		from functools import partial

		max_num_cores = os.cpu_count()

		if num_cpu > max_num_cores:
			raise ValueError('There are not enough cpus in the machine to run this action with ' + str(num_cpu) + ' cores.')
	print(text_color.GREEN + 'Inversion process has started..' + text_color.END)
	if num_cpu > 1:

		manager = multiprocessing.Manager()
		shared_results = manager.list()

		with multiprocessing.Pool(processes=num_cpu) as pool:

			process_item_partial = partial(_solv_cond_, cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
			lowerlimit=lower_limit_list , search_increment= search_start, acceptence_threshold = acceptence_threshold, results_list = shared_results, init_guess = None,
			transition_zone = transition_zone, water_solv=water_solv,comp_solv = comp_solv, melt_solv = melt_solv, comp_type = comp_type, comp_index = comp_index,
			low_value_threshold = low_value_threshold,sfd = simplify_fluid_density)

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
				water_solv=water_solv, comp_solv = comp_solv, melt_solv = melt_solv, comp_type = comp_type, comp_index = comp_index, low_value_threshold = low_value_threshold,
				sfd = simplify_fluid_density)
			
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

def _misfit(val, val_external, norm = 'log'):

	"""Internal function to calculate misfit in log10."""

	if norm == 'log':
		#misfit at log scale
		misf = np.log10(val) - np.log10(val_external)
	else:
		misf = val - val_external
	return misf

def _likelihood(val, val_external, sigma, norm = 'log'):

	"An internal function to calculate likelihood"

	misf = _misfit(val, val_external, norm = norm)
	misf_2 = -misf**2 / (2 * sigma**2)
	like = np.exp(misf_2)
	return like, misf_2
	

def _solv_MCMC_two_param(index, cond_list, object, initial_params, param_name_1, param_name_2, upper_limits,
	lower_limits, sigma_cond,proposal_stds,n_iter,burning,water_solv, comp_solv, continue_bool, vp_list = None, vs_list = None, sigma_vp = None, sigma_vs = None,
	adaptive_alg = True, ideal_acceptance_bounds = [0.2,0.3], adaptive_check_length = 1000, comp_index = [0,0], step_size_limits = None,
	transition_zone = False):
	
	"""
	MCMC external solver for the conductivity_metropolis_hastings_two_param function for parallelization purposes.
	Users should not call this function directly.
	"""

	if continue_bool[index] == True:

		#Using Metropolis-Hastings algorithm

		frac_bool = [False, False]

		if 'frac' in param_name_1:
			if param_name_1 != 'melt_fluid_mass_frac':
				frac_bool[0] = True
				comp_index_sub = 0

		if 'frac' in param_name_2:
			if param_name_2 != 'melt_fluid_mass_frac':
				frac_bool[1] = True
				comp_index_sub = 1

		if sum(frac_bool) == 2:

			raise KeyError('Currently only one of the parameters chosen can be modal compositional parameter.')
			
		if 'melt_fluid_mass_frac' in [param_name_1, param_name_2]:
		
			melt_solv = True
			
		else:
			
			melt_solv = False

		param_1_init, param_2_init = initial_params[index]
		(param_1_max,param_2_max) = upper_limits
		(param_1_min,param_2_min) = lower_limits

		current_params = np.array([param_1_init, param_2_init])

		#Initial setting of the parameters.
		if sum(frac_bool) == 1:

			if object.solid_phase_method == 2: #if mineral
				_comp_list = [object.quartz_frac[index], object.plag_frac[index], object.amp_frac[index], object.kfelds_frac[index], object.opx_frac[index], object.cpx_frac[index],
					object.mica_frac[index], object.garnet_frac[index], object.sulphide_frac[index], object.graphite_frac[index], object.ol_frac[index], object.sp_frac[index], object.rwd_wds_frac[index],
					object.perov_frac[index], object.mixture_frac[index], object.other_frac[index]]
			else: #if rock
				_comp_list = [object.granite_frac[index],object.granulite_frac[index],object.sandstone_frac[index],object.gneiss_frac[index],object.amphibolite_frac[index],
					object.basalt_frac[index],object.mud_frac[index],object.gabbro_frac[index],object.other_rock_frac[index]]

			comp_old = _comp_list[comp_index[comp_index_sub]]

			if frac_bool[0] == True:
				comp_list = _comp_adjust_(np.array(_comp_list), param_1_init, comp_old)
			else:
				comp_list = _comp_adjust_(np.array(_comp_list), param_2_init, comp_old)

			for idx_t in range(len(_comp_list)):

				if object.solid_phase_method == 2: #if mineral
					object.mineral_frac_list[idx_t][index] = comp_list[idx_t]
				else: #if rock
					object.rock_frac_list[idx_t][index] = comp_list[idx_t]

			water_solv = True

		#Executing the commands
		exec(f'object.{param_name_1}[{str(index)}] = param_1_init')
		exec(f'object.{param_name_2}[{str(index)}] = param_2_init')

		if water_solv == True:

			if transition_zone == False:
				object.mantle_water_distribute(method = 'index', sol_idx = index)
			else:
				object.transition_zone_water_distribute(method = 'index', sol_idx = index)
				
			if melt_solv == True:
				
				#to interpolation of fluid density so eos do not have to be solved at each iteration.
				try:
					water_index = [param_name_1,param_name_2].index('bulk_water')
					water_end = upper_limits[water_index] + (upper_limits[water_index] * 1000)
					water_end = water_end[0]
				except ValueError:
					water_end = 1e6

				object.calculate_density_fluid(sol_idx = index, method = 'array', interp_for_iter = True, water_start = 0, water_end = water_end)

		#Calculating the initial conductivity
		cond_init = object.calculate_conductivity(method = 'index', sol_idx = index)
		if (vp_list is not None) or (vs_list is not None):
			v_bulk_init, vp_init, vs_init = object.calculate_seismic_velocities(method = 'index',sol_idx = index)
			
		current_likelihood_cond, current_misf = _likelihood(cond_init, cond_list[index], sigma_cond[index])
		
		if vp_list is not None:
			current_likelihood_vp, misf_vp = _likelihood(vp_init, vp_list[index], sigma_vp[index])
		else:
			current_likelihood_vp = 1
			
		if vs_list is not None:
			current_likelihood_vs, misf_vs = _likelihood(vs_init, vs_list[index], sigma_vs[index])
		else:
			current_likelihood_vs = 1
			
		current_likelihood = current_likelihood_cond * current_likelihood_vp * current_likelihood_vs

		#empty arrays to fill it up with samples
		samples = []
		misfits_cond = []
		misfits_vp = []
		misfits_vs = []
		misfits_all_cond = []
		misfits_all_vp = []
		misfits_all_vs = []
		samples_all = []
		acceptance_rates = []
		accepted = 0
		print(text_color.GREEN + 'Monte-Carlo loop is started' + text_color.END)
		print(text_color.YELLOW + f'{n_iter} total samples.' + text_color.END)
		print(text_color.RED + f'{burning} burning samples.' + text_color.END)

		#loop for monte-carlo
		for _ in range(n_iter):
			
			#proposing the new parameters
			proposal = np.array(current_params)
			rand_node = int(np.floor(np.random.rand(1)*2.0)[0]) #randomnode generation 1 or 0
			randomgen = np.random.normal(0, proposal_stds[rand_node], size=1)[0] #proposal for randomwalk step
			proposal[rand_node] = current_params[rand_node] + randomgen
			
			continue_bounds = False
			if (proposal[0] > param_1_min[index]) and (proposal[0] < param_1_max[index]):
				if (proposal[1] > param_2_min[index]) and (proposal[1] < param_2_max[index]):
					continue_bounds = True
			else:
				proposal = current_params #if out of bounds go back to the previous parameters

			if continue_bounds == True:
				#adjusting for composition if needed...
				if comp_solv == True:

					if sum(frac_bool) == 1:

						if object.solid_phase_method == 2: #if mineral
							_comp_list = [object.quartz_frac[index], object.plag_frac[index], object.amp_frac[index], object.kfelds_frac[index], object.opx_frac[index], object.cpx_frac[index],
								object.mica_frac[index], object.garnet_frac[index], object.sulphide_frac[index], object.graphite_frac[index], object.ol_frac[index], object.sp_frac[index], object.rwd_wds_frac[index],
								object.perov_frac[index], object.mixture_frac[index], object.other_frac[index]]
						else: #if rock
							_comp_list = [object.granite_frac[index],object.granulite_frac[index],object.sandstone_frac[index],object.gneiss_frac[index],object.amphibolite_frac[index],
								object.basalt_frac[index],object.mud_frac[index],object.gabbro_frac[index],object.other_rock_frac[index]]

						comp_old = _comp_list[comp_index[comp_index_sub]]

						if frac_bool[0] == True:
							comp_list = _comp_adjust_(np.array(_comp_list), proposal[0], comp_old)
						else:
							comp_list = _comp_adjust_(np.array(_comp_list), proposal[1], comp_old)

						for idx_t in range(len(_comp_list)):

							if object.solid_phase_method == 2: #if mineral
								object.mineral_frac_list[idx_t][index] = comp_list[idx_t]
							else: #if rock
								object.rock_frac_list[idx_t][index] = comp_list[idx_t]

				#setting up the random parameter
				exec(f'object.{param_name_1}[{str(index)}] = proposal[0]')
				exec(f'object.{param_name_2}[{str(index)}] = proposal[1]')

				#distribute water if needed
				if water_solv == True:
					if transition_zone == False:
						object.mantle_water_distribute(method = 'index', sol_idx = index)
					else:
						object.transition_zone_water_distribute(method = 'index', sol_idx = index)

				proposed_cond = object.calculate_conductivity(method = 'index',sol_idx = index)
				if (vp_list is not None) or (vs_list is not None):
					v_bulk, proposed_vp, proposed_vs = object.calculate_seismic_velocities(method = 'index',sol_idx = index)
				
				
				proposed_likelihood_cond, misf_cond = _likelihood(proposed_cond, cond_list[index], sigma_cond[index])
				
				if vp_list is not None:
					proposed_likelihood_vp, misf_vp = _likelihood(proposed_vp, vp_list[index], sigma_vp[index], norm = 'log')
					print('####')
					print(proposed_vp, proposed_likelihood_vp,proposal[0],proposal[1],misf_vp)
					print(proposed_cond, proposed_likelihood_cond,proposal[0],proposal[1],misf_cond)
				else:
					proposed_likelihood_vp = 1
					misf_vp = 0
					
				if vs_list is not None:
					proposed_likelihood_vs, misf_vs = _likelihood(proposed_vs, vs_list[index], sigma_vs[index], norm = 'log')
				else:
					proposed_likelihood_vs = 1
					misf_vs = 0
				
				proposed_likelihood = proposed_likelihood_cond * proposed_likelihood_vs * proposed_likelihood_vp
				# print(proposed_likelihood_cond, proposed_vp, proposed_vs)
				# Calculate acceptance probability
				acceptance_ratio = proposed_likelihood / current_likelihood

				if np.random.rand() < acceptance_ratio:

					current_params = proposal
					current_likelihood = proposed_likelihood

					if _ > burning:
						samples.append(current_params)
						misfits_cond.append(misf_cond)
						misfits_vp.append(misf_vp)
						misfits_vs.append(misf_vs)
						accepted += 1

				if _ > burning:
					acceptance_rate = accepted / (_ - burning)
					acceptance_rates.append(acceptance_rate)
					misfits_all_cond.append(misf_cond)
					misfits_all_vp.append(misf_vp)
					misfits_all_vs.append(misf_vs)
					samples_all.append(current_params)
					if adaptive_alg == True:
						if (_ + 1) % adaptive_check_length == 0:
							if acceptance_rate <= ideal_acceptance_bounds[0]:
								proposal_stds[rand_node] = proposal_stds[rand_node] * 0.95
								if step_size_limits is not None:
									if proposal_stds[rand_node] > step_size_limits[rand_node]:
										proposal_stds[rand_node] = step_size_limits[rand_node]
								print(text_color.YELLOW + f'Step size (std) for random walk are decreased to {proposal_stds} - Acceptance Rate: {round(acceptance_rate,3)}- Completed :% {round((_/n_iter)*1e2)}' + text_color.END)
							elif acceptance_rate >= ideal_acceptance_bounds[1]:
								proposal_stds[rand_node] = proposal_stds[rand_node] * 1.05
								if step_size_limits is not None:
									if proposal_stds[rand_node] > step_size_limits[rand_node]:
										proposal_stds[rand_node] = step_size_limits[rand_node]
								print(text_color.RED + f'Step size (std) for random walk increased to {proposal_stds} - Acceptance Rate: {round(acceptance_rate,3)} - Completed :% {(_/n_iter)*1e2}' + text_color.END)
							else:
								print(text_color.GREEN + f'Acceptence rate is good size: - Acceptance Rate: {round(acceptance_rate,3)} - Completed :% {round((_/n_iter)*1e2)}' + text_color.END)
					else:
						if (_ + 1) % adaptive_check_length == 0:
							print(text_color.GREEN + f'Acceptance Rate: {round(acceptance_rate,3)}' + text_color.END)
	else:

		samples = np.array([None])
		acceptance_rates = np.array([None])
		misfits = np.array([None])
		samples_all = np.array([None])
		misfits_all = np.array([None])

		print(f'The value of {index}th index is below the dry conductivity of composition entered, therefore no solution is required.')
	
	misfits = [misfits_cond,misfits_vp,misfits_vs]
	misfits_all = [misfits_cond, misfits_vp, misfits_vs]
	
	return np.array(samples), np.array(acceptance_rates), misfits, np.array(samples_all), np.array(misfits_all)

def conductivity_metropolis_hastings_two_param(object, cond_list, initial_params, param_name_1, param_name_2, upper_limits,
	lower_limits, sigma_cond, proposal_stds,n_iter, vp_list = None, vs_list = None, sigma_vs= None, sigma_vp = None, burning = 0, 
	transition_zone = False, num_cpu = 1, **kwargs):

	"""
	Perform Metropolis-Hastings MCMC inversion for electrical conductivity using two model parameters.

	This function uses a stochastic sampling approach to estimate the posterior distribution of two 
	input parameters based on observed conductivity data. It is useful for exploring uncertainty 
	and trade-offs between parameters in conductivity models.

	Parameters
	----------
	object : object
		A pide model instance for calculating conductivity.
	cond_list : array-like
		Observed conductivity values to fit [S/m].
	initial_params : list or array-like
		Initial values for the two parameters to invert.
	param_name_1 : str
		Name of the first parameter to invert (must be an attribute of `object`).
	param_name_2 : str
		Name of the second parameter to invert (must be an attribute of `object`).
	upper_limits : tuple of floats
		Upper bounds for the parameter search space of first and second param.
	lower_limits : tuple of floats
		Lower bounds for the parameter search space of first and second param.
	sigma_cond : float or array-like
		Standard deviation or error for conductivity observations in logarithm of conductivity [S/m].
		Do not enter this as a percentage error, you need to enter as absolute error in log-conductivity.
	proposal_stds : list or array-like
		Standard deviations for proposal distribution of the two parameters.
		This acts as a initial step size in for the random search.
	n_iter : int
		Number of iterations for the MCMC chain.
	burning : int, optional
		Number of initial iterations to discard as burn-in (default is 0).
	transition_zone : bool, optional
		If True, use transition zone water distribution functions (default is False).
	num_cpu : int, optional
		Number of CPU cores to use for parallel computation (default is 1).
	save_distr : bool, optional
		If True, saves the MCMC samples to disk (default is False).
	distr_file_names : str, optional
		Base name for saved distribution files (default is 'distribution_solution').
	adaptive_alg : bool, optional
		If True, enables adaptive adjustment of proposal standard deviations based on acceptance rate.
	ideal_acceptance_bounds : list of float, optional
		Target acceptance rate bounds for adaptive algorithm (default is [0.2, 0.3]).
	adaptive_check_length : int, optional
		Number of iterations between check for adaptive algorithm (default is 1000).
	step_size_limits : list of float, optional
		Minimum and maximum bounds for the adaptive proposal step sizes.

	Returns
	-------
	Sample distribution (accepted): array-like
		Array of sampled accepted parameter values of shape (n_accepted_samples, 2).
	acceptance_rate : array-like
		Acceptenca rate record over the sampling.
	misfits : array-like
		misfits of the accepted distribution.
	Sample distribution (all) : array-like
		Array of sampled all parameter values of shape (n_iter - burning, 2).
	misfits_all: array-like
		
	sample_distr, acceptance_rates, misfits, samples_all, misfits_all

	Examples
	--------
	samples, acceptance_rates, misfits, samples_all, misfits_all = conductivity_metropolis_hastings_two_param(object = p_obj, cond_list = [0.1,0.1],
	initial_params = [[200,0.25]],param_name_1 = 'bulk_water',
	param_name_2= "melt_fluid_mass_frac", upper_limits = (2000,0.5),
	lower_limits = (0,0), sigma_cond = [0.1,0.1],proposal_stds=[200,0.25]
	,n_iter = 2e5, burning = 1e4, transition_zone = False,num_cpu = 1,adaptive_alg = True,
	step_size_limits = [25000,0.5])
	
	"""

	#Pre-checks for if
	if object.solid_phase_method == 2:
		object.set_mineral_water(ol = 0, opx = 0, cpx = 0, garnet = 0, mica = 0, amp = 0,
		quartz = 0, plag = 0, kfelds = 0, sulphide = 0, graphite = 0, sp = 0, rwd_wds = 0,
		perov = 0, mixture = 0, other = 0)
	elif object.solid_phase_method == 1:
		object.set_rock_water(granite = 0, granulite = 0, sandstone = 0, gneiss = 0,
		amphibolite = 0, basalt = 0, mud = 0, gabbro = 0, other_rock = 0)

	cond_check = object.calculate_conductivity()

	save_distr = kwargs.pop('save_distr',False)
	distr_file_names = kwargs.pop('distr_file_names','distribution_solution')
	adaptive_alg = kwargs.pop('adaptive_alg', True)
	ideal_acceptance_bounds = kwargs.pop('ideal_acceptance_bounds',[0.2,0.3])
	adaptive_check_length = kwargs.pop('adaptive_check_length', 1000)
	step_size_limits = kwargs.pop('step_size_limits',None)

	#Pre checks for the input parameters.
	if type(ideal_acceptance_bounds) == list:
		if len(ideal_acceptance_bounds) == 2:
			pass
		else:
			raise ValueError(f'ideal_acceptance_bounds has to be a list containing two values. Currently it is {ideal_acceptance_bounds}')
	else:
		raise ValueError(f'ideal_acceptance_bounds has to be a list containing two values. Currently it is {ideal_acceptance_bounds}')

	try:
		exec(f'object.{param_name_1}')
		exec(f'object.{param_name_2}')
	except AttributeError:
		raise AttributeError(f'There is no such parameter name {param_name_1} or {param_name_2} for the pide object.')

	continue_bool = []
	for i in range(len(cond_list)):
		continue_bool.append(cond_list[i] > cond_check[i])

	if burning >= n_iter:

		raise ValueError('Burning samples cannot be larger than the total iteration number (n_iter).')

	if len(cond_list) == len(initial_params) == len(upper_limits[0]) == len(lower_limits[0]) == len(sigma_cond) == len(object.T):
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

	comp_solv = False
	water_solv = False
	comp_type = None
	comp_index = []
	comp_type_list = []

	if any('water' in xx for xx in param_names) == True:

		if 'bulk_water' in param_names:
			
			water_solv = True
			#setting the object.bulk_water as same length as T if that has not done already...
			if len(getattr(object,param_names[param_names.index('bulk_water')])) != len(object.T):
				object.set_bulk_water(0.0)
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')

	if any('melt' in xx for xx in param_names) == True:

		water_solv = True
		
		for ii in range(2):
			if len(getattr(object,param_names[ii])) != len(object.T):
				object.set_parameter(param_names[ii], 0.0)

	if any('frac' in xx for xx in param_names) == True:

		comp_solv = True
		water_solv = True

		for ii in range(len(param_names)):

			if param_names[ii] != 'melt_fluid_mass_frac':

				if param_names[ii] in min_list:

					comp_type = 'mineral'
					comp_type_list.append(comp_type)
					comp_index.append(min_list.index(param_names[ii]))

				elif param_names[ii] in rock_list:

					comp_type = 'rock'
					comp_type_list.append(comp_type)
					comp_index.append(rock_list.index(param_names[ii]))

				else:

					comp_type = None
					comp_type_list.append(comp_type)
					comp_index.append(None)

			else:

				comp_type = None
				comp_type_list.append(comp_type)
				comp_index.append(None)

			if comp_type is not None:
				if len(getattr(object,param_names[ii])) != len(object.T):
					object.set_parameter(param_names[ii], 0.0)

		if (('mineral' in comp_type_list) == True) and (('rock' in comp_type_list) == True):
			raise ValueError('The user cannot enter both rock and mineral as the inversion parameter. Choose only one.')

	#The last check for setting up other parameters.
	for ii in range(len(param_names)):
		if len(getattr(object,param_names[ii])) != len(object.T):
			object.set_parameter(param_names[ii], 0.0)

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
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond, proposal_stds = proposal_stds , n_iter= n_iter, burning = burning, vp_list = vp_list, vs_list = vs_list,
			sigma_vp = sigma_vp, sigma_vs = sigma_vs,
			water_solv = water_solv, comp_solv = comp_solv, comp_index = comp_index, continue_bool = continue_bool, adaptive_alg = adaptive_alg, adaptive_check_length = adaptive_check_length,
			step_size_limits = step_size_limits, ideal_acceptance_bounds = ideal_acceptance_bounds)

			c = pool.map(process_item_partial, index_list)

		sample_distr = [x[0] for x in c]
		acceptance_rates = [x[1] for x in c]
		misfits = [x[2] for x in c]
		samples_all = [x[3] for x in c]
		misfits_all = [x[4] for x in c]

	else:

		sample_distr = []
		acceptance_rates = []
		misfits = []
		samples_all = []
		misfits_all = []

		for idx in range(0,len(index_list)):
			
			c = _solv_MCMC_two_param(index = index_list[idx], object = object, cond_list = cond_list, initial_params = initial_params, param_name_1 = param_name_1, param_name_2= param_name_2,
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond, proposal_stds = proposal_stds , n_iter= n_iter, burning = burning, vp_list = vp_list, vs_list = vs_list,
			sigma_vp = sigma_vp, sigma_vs = sigma_vs,
			water_solv = water_solv, comp_solv = comp_solv, comp_index = comp_index, continue_bool = continue_bool, adaptive_alg = adaptive_alg, adaptive_check_length = adaptive_check_length,
			step_size_limits = step_size_limits, ideal_acceptance_bounds = ideal_acceptance_bounds)
			
			sample_distr.append(c[0])
			acceptance_rates.append(c[1])
			misfits.append(c[2])
			samples_all.append(c[3])
			misfits_all.append(c[4])

	if save_distr == True:

		from pide.utils.utils import save_h5_files

		array_names_idx = list(range(len(sample_distr)))
		array_names_idx = [str(element) for element in array_names_idx]

		save_h5_files(array_list=sample_distr, array_names=array_names_idx, file_name=distr_file_names + '_distr.h5')
		save_h5_files(array_list=acceptance_rates, array_names=array_names_idx, file_name=distr_file_names + '_acceptance.h5')
		save_h5_files(array_list=misfits, array_names=array_names_idx, file_name=distr_file_names + '_misfit.h5')
		save_h5_files(array_list=samples_all, array_names=array_names_idx, file_name=distr_file_names + '_distr_all.h5')
		save_h5_files(array_list=misfits_all, array_names=array_names_idx, file_name=distr_file_names + '_misfit_all.h5')

	return sample_distr, acceptance_rates, misfits, samples_all, misfits_all