#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from pide.geodyn.geotherm import calculate_hasterok2011_geotherm, T_Katsura_2022_Adiabat
from .utils.utils import text_color, check_type
import copy
from pide.utils.utils import save_h5_files

def _comp_adjust_(_comp_list, comp_alien, comp_old,final = False):

	"""A method to adjust composition of one mineral/rock without considering the replacement weights
	"""

	if final == False:
		ratio = (comp_alien - comp_old) / (np.sum(_comp_list) - comp_old)
	else:
		ratio = (comp_alien - comp_old) / (np.sum(_comp_list,axis = 0) - comp_old)

	comp_list = _comp_list - (_comp_list * ratio)

	return comp_list

def _solv_cond_(index, cond_list, object, param, upperlimit, lowerlimit, acceptence_threshold, results_list = None,
	init_guess = None, transition_zone = False, water_solv=False, comp_solv=False, melt_solv=False, comp_type = None, comp_index = None, low_value_threshold = None,
	sfd = False,init_guess_preiter = True):

	"""Simple line-search solver to solve conductivities. This function should not be called directly.
	"""
	
	if init_guess_preiter == True:
		if results_list is not None:
			if len(results_list) == 0:
				init_guess = None
			else:
				init_guess = results_list[-1]
				
	else:
		init_guess = None
	
	search_increment = (upperlimit[index] - lowerlimit[index]) / 8.0
	param_search_array = np.arange(lowerlimit[index], upperlimit[index] , search_increment)
		
	restart = True
	init_search_increment = search_increment.copy()
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
						sol_param = lowerlimit[index]
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
						sol_param = upperlimit[index] #equivalent to upper limit

						restart = False
						break
					else:
						pass
					
	if results_list is not None:
		results_list.append(sol_param)

	return sol_param, residual


def conductivity_solver_single_param(object, cond_list, param_name,
	upper_limit_list, lower_limit_list, acceptence_threshold, cond_err = None, transition_zone = False, simplify_fluid_density = False,
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
	acceptence_threshold = 1, cond_err = None, transition_zone = False, simplify_fluid_density = False,
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
	init_guess_preiter = kwargs.pop('init_guess_preiter', True)

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
			lowerlimit=lower_limit_list , acceptence_threshold = acceptence_threshold, results_list = shared_results, init_guess = None,
			transition_zone = transition_zone, water_solv=water_solv,comp_solv = comp_solv, melt_solv = melt_solv, comp_type = comp_type, comp_index = comp_index,
			low_value_threshold = low_value_threshold,sfd = simplify_fluid_density,init_guess_preiter = init_guess_preiter)

			c = pool.map(process_item_partial, index_list)

		c_list = [x[0] for x in c]
		residual_list= [x[1] for x in c]

		c_list = np.array(c_list)
		residual_list = np.array(residual_list)

	else:

		c_list = np.zeros(len(index_list))
		residual_list = np.zeros(len(index_list))
		
		for idx in range(0,len(index_list)):

			if idx > 0:
				if init_guess_preiter == False:
					init_guess_ = None
				else:
					init_guess_ = c_list[idx-1]
			else:
				init_guess_ = None

			c = _solv_cond_(index = index_list[idx], cond_list = cond_list, object = object, param = param_name, upperlimit = upper_limit_list,
				lowerlimit=lower_limit_list , acceptence_threshold = acceptence_threshold, results_list= None, init_guess = init_guess_, transition_zone = transition_zone,
				water_solv=water_solv, comp_solv = comp_solv, melt_solv = melt_solv, comp_type = comp_type, comp_index = comp_index, low_value_threshold = low_value_threshold,
				sfd = simplify_fluid_density, init_guess_preiter = init_guess_preiter)
			
			c_list[idx] = c[0]
			residual_list[idx] = c[1]
			
	print(text_color.CYAN + 'Inversion process has ended..' + text_color.END)
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
			
		elif sum(frac_bool) > 1:
		
			raise ValueError('The user cannot change more than 1 modal fraction value.')

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
					water_end = upper_limits[water_index]
					water_end = water_end[0]
				except ValueError:
					water_end = 1e5

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
					proposed_likelihood_vp, misf_vp = _likelihood(proposed_vp, vp_list[index], sigma_vp[index], norm = 'linear')
				else:
					proposed_likelihood_vp = 1
					misf_vp = 0
					
				if vs_list is not None:
					proposed_likelihood_vs, misf_vs = _likelihood(proposed_vs, vs_list[index], sigma_vs[index], norm = 'linear')
				else:
					proposed_likelihood_vs = 1
					misf_vs = 0
				
				proposed_likelihood = proposed_likelihood_cond * proposed_likelihood_vs * proposed_likelihood_vp
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

def metropolis_hastings_two_param(object, cond_list, initial_params, param_name_1, param_name_2, upper_limits,
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

def _fractions_to_unconstrained(f_pyx, f_lherz, eps=1e-6):
	"""Real (f_pyx, f_lherz) -> two unconstrained reals that can take any
	value. f_dun = 1 - f_pyx - f_lherz is implied."""
	f_dun = 1.0 - f_pyx - f_lherz
	u1 = np.clip(f_dun, eps, 1 - eps)
	remaining = f_pyx + f_lherz
	u2 = np.clip(f_pyx / remaining, eps, 1 - eps) if remaining > eps else 0.5
	return np.log(u1 / (1 - u1)), np.log(u2 / (1 - u2))


def _unconstrained_to_fractions(x1, x2):
	"""Two unconstrained reals -> valid (f_pyx, f_lherz). ANY input gives a
	pair satisfying f_pyx >= 0, f_lherz >= 0, f_pyx + f_lherz <= 1, so the
	triangle constraint can never be violated and never needs checking."""
	u1 = 1.0 / (1.0 + np.exp(-x1))
	u2 = 1.0 / (1.0 + np.exp(-x2))
	remaining = 1.0 - u1
	return remaining * u2, remaining * (1.0 - u2)
	
def _solv_MCMC_column(index, object, depths, moho_depth,
	cond_list, vp_list, vs_list,
	sigma_cond, sigma_vp, sigma_vs,
	initial_SHF, initial_lab_temp,
	initial_params, param_names,
	upper_limits, lower_limits,
	SHF_bounds,
	proposal_stds, n_iter, burning, lab_depth = None, sigma_lab = None,
	melt_thermodyn=False, melt_thermodyn_interp=None,
	adaptive_alg=True, ideal_acceptance_bounds=[0.2, 0.3],
	adaptive_check_length=1000, step_size_limits=None,
	param_priors=None, SHF_prior=None, lab_temp = 1350,
	invert_lab_temp = False, lab_temp_bounds = None,
	composition = None, max_widen_attempts = 3,
	record_names = None,
	**kwargs):
	"""
	MCMC column solver for geotherm-based inversion.

	Samples SHF as a scalar parameter, plus depth-varying parameters
	(e.g., bulk_water, bulk_xfe) along a depth column. Temperature is
	derived from the geotherm function using SHF and a LAB temperature.

	If invert_lab_temp is False (default), lab_temp is held FIXED at the
	value given by the `lab_temp` argument (Celsius), exactly as before —
	this is a backward-compatible default, nothing else changes.

	If invert_lab_temp is True, LAB temperature becomes a SECOND sampled
	scalar parameter alongside SHF, using `initial_lab_temp` as its
	starting value (per column, same shape as initial_SHF) and
	`lab_temp_bounds` as its (min, max) bounds (same shape as SHF_bounds).
	The fixed `lab_temp` argument is ignored in this case.

	IMPORTANT: proposal_stds ordering changes depending on invert_lab_temp:
		invert_lab_temp=False: [SHF_std, param0_std, param1_std, ...]
		invert_lab_temp=True:  [SHF_std, lab_temp_std, param0_std, param1_std, ...]

	record_names lets you additionally record ANY other object attribute
	(e.g. ['ol_frac', 'ol_xfe', 'garnet_frac']) at every depth, alongside
	the sampled parameters, without those attributes being sampled
	themselves — they're just read off the object and stored, wherever
	they end up after everything else that iteration set. The recorded
	dicts are always the LAST two items in the return tuple, regardless
	of which other flags (invert_lab_temp, melt_thermodyn) are set, so
	results[-2] / results[-1] always work.

	Parameters
	----------
	index : int
		Column index for parallel dispatch.
	object : pide object
		Pre-configured pide instance with length = len(depths).
	depths : array
		Depth nodes in km.
	moho_depth : float
		Moho depth in km (fixed).
	cond_list : array
		Observed conductivity at each depth node [S/m].
	vp_list : array or None
		Observed Vp at each depth node [km/s].
	vs_list : array or None
		Observed Vs at each depth node [km/s].
	lab_depth: float or None
		Observed Lab depth [km] from independent information. This is used
		to fit to the estimations.
	sigma_cond : array
		Conductivity uncertainty at each depth (log space).
	sigma_vp : array or None
		Vp uncertainty at each depth [km/s].
	sigma_vs : array or None
		Vs uncertainty at each depth [km/s].
	sigma_lab : float or None
		Lab depth uncertainty [km].
	initial_SHF : array
		Initial surface heat flow [mW/m^2], per column.
	initial_lab_temp : array
		Initial LAB temperature [Celsius], per column. Only used if
		invert_lab_temp=True; otherwise ignored (pass anything, e.g. None).
	initial_params : array
		Shape (n_depths, n_params). Initial values for depth-varying params.
	param_names : list of str
		Names of depth-varying parameters (e.g., ['bulk_water', 'bulk_xfe']).
	upper_limits : tuple of arrays
		Upper bounds per parameter. Each array length = n_depths.
	lower_limits : tuple of arrays
		Lower bounds per parameter. Each array length = n_depths.
	SHF_bounds : tuple
		(min_SHF, max_SHF), per column.
	proposal_stds : list
		Step sizes. See ordering note above.
	n_iter : int
		Total MCMC iterations.
	burning : int
		Burn-in iterations.
	melt_thermodyn : bool
		If True, calculate melt from Katz lookup.
	melt_thermodyn_interp : RegularGridInterpolator or None
		3D interpolator (T_celsius, water_wt%, P_GPa) -> melt fraction.
	adaptive_alg : bool
		Enable adaptive step sizes.
	ideal_acceptance_bounds : list
		Target acceptance rate [low, high].
	adaptive_check_length : int
		Iterations between adaptive checks.
	step_size_limits : list or None
		[min, max] step size per parameter, same ordering as proposal_stds.
	param_priors : list or None
		Per depth-varying parameter. Each is None or (mean_array, sigma_array).
	SHF_prior : tuple or None
		(mean, sigma) for Gaussian SHF prior.
	lab_temp : float
		Fixed LAB temperature in Celsius, used when invert_lab_temp=False.
		Default 1350.
	invert_lab_temp : bool
		If True, sample LAB temperature as a free scalar parameter instead
		of holding it fixed. Default False.
	lab_temp_bounds : tuple or None
		(min_lab_temp, max_lab_temp), per column. Required if
		invert_lab_temp=True; ignored otherwise.
	record_names : list of str or None
		Names of additional object attributes to record at every depth
		(e.g. ['ol_frac', 'ol_xfe']), read off the object as-is, not
		sampled. Default None (records nothing extra).
	"""

	melt_frac_limit = kwargs.pop('melt_frac_limit',0.005)
	comp_index = kwargs.pop('comp_index',[0] * len(param_names))

	#deep copy object to not confuse multiprocessing workers.
	object = copy.deepcopy(object)
	
	widen_count = 0

	#determining length of the parametrisation
	n_depths = len(depths)
	n_params = len(param_names)

	frac_bool = [False] * n_params
	comp_index_sub = None

	for ii in range(n_params):
		if 'frac' in param_names[ii]:
			if param_names[ii] != 'melt_fluid_mass_frac':
				frac_bool[ii] = True
				comp_index_sub = ii

	# Temperature and pressure are controlled by the geotherm, not directly sampled
	for name in param_names:
		if name in ['T', 'p']:
			raise ValueError(f"'{name}' cannot be a free parameter in column mode. "
				f"Temperature is controlled by SHF and pressure is derived from the geotherm.")

	if ('f_pyx' in param_names) or ('f_lherz' in param_names):
		triangle_calculation = True
		if 'bulk_xfe' in param_names:
			raise AttributeError('Both triangle composition (f_pyx or f_lherz) and bulk_xfe can not be considered to vary with composition at the same time.')
	else:
		triangle_calculation = False

	if invert_lab_temp == True and lab_temp_bounds is None:
		raise ValueError('lab_temp_bounds must be provided when invert_lab_temp=True.')

	if record_names is None:
		record_names = []

	# --- Scalars: SHF is always sampled; lab_temp is optionally a second
	# sampled scalar. Each scalar gets n_params "slots" of sampling weight,
	# same relative emphasis SHF always had -- this generalizes the old
	# single-scalar setup without changing behavior when invert_lab_temp=False.
	scalar_names = ['SHF']
	scalar_slots = [n_params]
	scalar_bounds_list = [SHF_bounds]
	current_scalars = [initial_SHF[index]]

	if invert_lab_temp == True:
		scalar_names.append('lab_temp')
		scalar_slots.append(n_params)
		scalar_bounds_list.append(lab_temp_bounds)
		current_scalars.append(initial_lab_temp[index])

	current_scalars = np.array(current_scalars, dtype=float)
	n_scalars = len(scalar_names)
	n_scalar_total = sum(scalar_slots)
	scalar_slot_bounds = np.cumsum([0] + scalar_slots)

	n_param_total = n_params * n_depths
	n_joint_total = n_depths if triangle_calculation == True else 0
	n_total = n_param_total + n_scalar_total + n_joint_total

	def _get_scalar_idx(dim):
		return int(np.searchsorted(scalar_slot_bounds, dim, side='right') - 1)

	def _get_step_idx(dim):
		if dim < n_scalar_total:
			return _get_scalar_idx(dim)
		elif dim < n_scalar_total + n_param_total:
			dim_offset = dim - n_scalar_total
			return n_scalars + dim_offset // n_depths
		else:
			return n_scalars + n_params  # dedicated step index for the joint move

	# Deep copy mutable inputs
	proposal_stds = list(proposal_stds)
	#adding another step for the joint step if triangle_calculation is True.
	if triangle_calculation == True:
		proposal_stds = list(proposal_stds) + [np.mean(proposal_stds[-2:])]
	if param_priors is not None:
		param_priors = copy.deepcopy(param_priors)

	current_depth_params = np.array(initial_params, dtype=float)  # shape (n_depths, n_params)
	if current_depth_params.ndim == 1:
		current_depth_params = current_depth_params.reshape(-1, 1)  # (n_depths, 1)

	initial_params = np.array(initial_params, dtype=float)
	if initial_params.ndim == 1:
		initial_params = initial_params.reshape(-1, 1)

	# --- Determine which params need water distribution ---
	water_solv = 'bulk_water' in param_names

	melt_solv = 'melt_fluid_mass_frac' in param_names
	if melt_thermodyn is True:
		melt_solv = True

	# --- Calculate initial geotherm and set on object ---
	def _update_geotherm(shf_val, lab_temp_val):
		"""Calculate geotherm and interpolate onto depth nodes."""
		T_full, depth_full, p_full, idx_LAB = calculate_hasterok2011_geotherm(
			SHF=shf_val, T_0=0, max_depth=depths[-1] + 10,
			moho=moho_depth[index], adiabat=True,
			thermal_lab = True, thermal_lab_temp = lab_temp_val)

		# Interpolate onto our depth nodes
		T_interp_func = interp1d(depth_full, T_full, kind='linear', fill_value='extrapolate')
		P_interp_func = interp1d(depth_full, p_full, kind='linear', fill_value='extrapolate')

		T_at_depths = T_interp_func(depths)
		P_at_depths = P_interp_func(depths)
		LAB = depth_full[idx_LAB]

		return T_at_depths, P_at_depths, LAB

	#Generating the initial temperature from initial SHF (and lab_temp if inverted) distribution.
	_lab_temp_for_geotherm = current_scalars[1] if invert_lab_temp else lab_temp
	T_init, P_init, LAB = _update_geotherm(shf_val=current_scalars[0], lab_temp_val=_lab_temp_for_geotherm)
	object.set_temperature(T_init)
	object.set_pressure(P_init)

	#Setting the composition if defined, recasting into temperature and pressure length if single.
	if composition is not None:
		object.set_composition_solid_mineral(**composition)
	else:
		object.set_composition_solid_mineral(reval = True)

	#Setting all the parameters defined...
	for ii in range(n_params):
		getattr(object, param_names[ii])[:n_depths] = current_depth_params[:, ii]

	# converting the f_pyx/f_lherz columns into unconstrained coordinates ---
	# From here on, these two columns of current_depth_params hold logit-space
	# values (any real number), NOT fractions. They get converted back to real
	# fractions only where the physics needs them, and before samples are stored.
	if triangle_calculation == True:
		idx_fp = param_names.index('f_pyx')
		idx_fl = param_names.index('f_lherz')
		for iz in range(n_depths):
			x1, x2 = _fractions_to_unconstrained(
				current_depth_params[iz, idx_fp], current_depth_params[iz, idx_fl])
			current_depth_params[iz, idx_fp] = x1
			current_depth_params[iz, idx_fl] = x2

	#if bulk xFe is changed distributing iron among defined minerals.
	if 'bulk_xfe' in param_names:
		if triangle_calculation == False:
			object.mantle_xfe_distribute()
		else:
			raise AttributeError('xFe and triangle gibbs calculation method cannot be used at the same time.')
	else:
		if triangle_calculation == True:
			object.calculate_composition_modulii_from_triangle(method = 'array')

	#Calculating thermodynamic calculation of partial melting is True.
	current_melt = None

	if melt_thermodyn and melt_thermodyn_interp is not None:

		idx_water = param_names.index('bulk_water') if 'bulk_water' in param_names else None
		current_melt = np.zeros(n_depths)
		for iz in range(n_depths):
			T_C = object.T[iz] - 273.15
			if idx_water is not None:
				water_wt = current_depth_params[iz, idx_water] * 1e-4
			else:
				water_wt = object.bulk_water[iz] * 1e-4
			melt_val = float(melt_thermodyn_interp([T_C, water_wt, object.p[iz]]))
			if melt_val < melt_frac_limit:
				melt_val = 0.0
			current_melt[iz] = melt_val
		object.melt_fluid_mass_frac[:n_depths] = current_melt

	# --- Set up fluid density interpolation for melt calculations ---
	if water_solv == True:

		object.mantle_water_distribute()

		if (melt_solv == True) or (melt_thermodyn == True):

			#to interpolation of fluid density so eos do not have to be solved at each iteration.
			try:
				idx_water = param_names.index('bulk_water')
				water_end = upper_limits[idx_water][index]
				water_end = np.amax(water_end)
			except ValueError:
				water_end = 1e5

			object.calculate_density_fluid(sol_idx = index, method = 'array',
			interp_for_iter = True, water_start = 0, water_end = water_end)


	if (vp_list is not None) or (vs_list is not None):
		v_bulk_init, vp_init, vs_init = object.calculate_seismic_velocities(method = 'array')

	#Calculating the initial conductivity
	if cond_list is not None:
		cond_init = object.calculate_conductivity(method = 'array')

	if cond_list is not None:
		current_likelihood_cond, misf_cond = _likelihood(cond_init, cond_list[index], sigma_cond[index])
		current_likelihood_cond = np.sum(current_likelihood_cond)
	else:
		current_likelihood_cond = 1
		misf_cond = np.zeros(len(object.T))

	if vp_list is not None:
		current_likelihood_vp, misf_vp = _likelihood(vp_init, vp_list[index], sigma_vp[index], norm = 'linear')
		current_likelihood_vp = np.sum(current_likelihood_vp)
	else:
		current_likelihood_vp = 1
		misf_vp = np.zeros(len(object.T))
	if vs_list is not None:
		current_likelihood_vs, misf_vs = _likelihood(vs_init, vs_list[index], sigma_vs[index], norm = 'linear')
		current_likelihood_vs = np.sum(current_likelihood_vs)
	else:
		current_likelihood_vs = 1
		misf_vs = np.zeros(len(object.T))
		
	if lab_depth is not None:
		current_likelihood_lab, misf_lab = _likelihood(LAB, lab_depth, sigma_lab, norm = 'linear')
	else:
		current_likelihood_lab = 1
		misf_lab = 0.0

	current_prior_log = 0.0

	#SHF prior:
	if SHF_prior is not None:
		current_prior_log += -0.5 * ((current_scalars[0] - SHF_prior[index][0]) / SHF_prior[index][1])**2

	#DO THIS LATER, does not work now.
	if param_priors is not None:
		for ii in range(n_params):
			if param_priors[ii] is not None:
				prior_mean = param_priors[ii][0]   # array length n_depths
				prior_sigma = param_priors[ii][1]   # array length n_depths

				current_prior_log += np.sum(-0.5 * ((initial_params[:, ii] - prior_mean) / prior_sigma)**2)

	current_likelihood = np.exp(np.sum(misf_cond) + np.sum(misf_vp) + np.sum(misf_vs) + misf_lab + current_prior_log)

	samples_SHF = []
	samples_lab_temp = []
	samples_lab = []
	samples_temp = []
	samples_depth_params = []
	misfits_cond = []
	misfits_vp = []
	misfits_vs = []
	misfits_lab = []
	misfits_all_cond = []
	misfits_all_vp = []
	misfits_all_vs = []
	misfits_all_lab = []
	samples_SHF_all = []
	samples_lab_temp_all = []
	samples_lab_all = []
	samples_temp_all = []
	samples_depth_params_all = []
	acceptance_rates = []
	melt_samples = []
	melt_samples_all = []
	
	accepted = 0

	# --- Additional recorded object attributes (not sampled, just read off
	# the object and stored alongside everything else). ---
	samples_record = {name: [] for name in record_names}
	samples_record_all = {name: [] for name in record_names}

	# --- Bounds ---
	param_mins_depth = np.array([lower_limits[i][index] for i in range(n_params)])  # (n_params, n_depths)
	param_maxs_depth = np.array([upper_limits[i][index] for i in range(n_params)])

	# Ensure 2D even with single parameter
	if param_mins_depth.ndim == 1:
		param_mins_depth = param_mins_depth.reshape(1, -1)
		param_maxs_depth = param_maxs_depth.reshape(1, -1)

	if triangle_calculation == True:
		# these columns are no longer fractions, so their user-supplied
		# [0,1] bounds no longer apply -- replace with a wide unconstrained
		# range (±6 in logit space covers fractions ~0.0025 to ~0.9975)
		f_pyx_real_floor = param_mins_depth[idx_fp, :].copy()    # your actual 0.05, saved
		f_lherz_real_floor = param_mins_depth[idx_fl, :].copy()
		param_mins_depth[idx_fp, :] = -6.0
		param_maxs_depth[idx_fp, :] = 6.0
		param_mins_depth[idx_fl, :] = -6.0
		param_maxs_depth[idx_fl, :] = 6.0

	def _to_real_fractions(params):
		"""Convert the stored logit-space f_pyx/f_lherz columns back to real
		fractions, for saving samples. Every other column (e.g. bulk_water)
		is copied through untouched, and nothing is converted at all if the
		triangle method isn't in use."""
		out = params.copy()
		if triangle_calculation == True:
			for iz in range(n_depths):
				fp, fl = _unconstrained_to_fractions(out[iz, idx_fp], out[iz, idx_fl])
				out[iz, idx_fp] = fp
				out[iz, idx_fl] = fl
		return out

	print(text_color.GREEN + 'Monte-Carlo loop is started' + text_color.END)
	print(text_color.YELLOW + f'{n_iter*2} total samples, {n_iter} minimum samples.' + text_color.END)
	print(text_color.RED + f'{burning} burning samples.' + text_color.END)
	if invert_lab_temp == True:
		print(text_color.YELLOW + 'LAB temperature is being inverted as a free scalar parameter.' + text_color.END)
	if len(record_names) > 0:
		print(text_color.YELLOW + f'Additionally recording: {record_names}' + text_color.END)

	n_reject_bounds = 0
	n_reject_nan = 0
	n_reject_likelihood = 0

	n_step_dims = n_scalars + n_params + (1 if triangle_calculation == True else 0)
	n_attempted_per_dim = [0] * n_step_dims
	n_accepted_per_dim = [0] * n_step_dims

	for _ in range(n_iter*2):

		# Pick random dimension
		
		rand_dim = np.random.randint(n_total)
		step_idx = _get_step_idx(rand_dim)
		n_attempted_per_dim[step_idx] += 1

		randomgen = np.random.normal(0, proposal_stds[step_idx])

		# Copy current state
		proposed_scalars = current_scalars.copy()
		proposed_depth_params = current_depth_params.copy()
		continue_bounds = True

		if rand_dim < n_scalar_total:
			s_idx = _get_scalar_idx(rand_dim)
			proposed_scalars[s_idx] = current_scalars[s_idx] + randomgen
			lo, hi = scalar_bounds_list[s_idx][index]
			if proposed_scalars[s_idx] < lo or proposed_scalars[s_idx] > hi:
				continue_bounds = False

		elif rand_dim < n_scalar_total + n_param_total:
			dim_offset = rand_dim - n_scalar_total
			param_idx = dim_offset // n_depths
			depth_idx = dim_offset % n_depths

			proposed_depth_params[depth_idx, param_idx] += randomgen

			if (proposed_depth_params[depth_idx, param_idx] < param_mins_depth[param_idx, depth_idx] or
				proposed_depth_params[depth_idx, param_idx] > param_maxs_depth[param_idx, depth_idx]):
				continue_bounds = False

			if triangle_calculation == True and continue_bounds == True and param_names[param_idx] in ('f_pyx', 'f_lherz'):
				fp_check, fl_check = _unconstrained_to_fractions(
					proposed_depth_params[depth_idx, idx_fp],
					proposed_depth_params[depth_idx, idx_fl])
				if fp_check < f_pyx_real_floor[depth_idx] or fl_check < f_lherz_real_floor[depth_idx]:
					continue_bounds = False

		else:
			# joint f_pyx/f_lherz move: step BOTH logit values together in
			# one iteration, fixing f_lherz's structural disadvantage in
			# the stick-breaking transform (see earlier discussion).
			joint_depth_idx = rand_dim - n_scalar_total - n_param_total

			eps_fp = np.random.normal(0, proposal_stds[step_idx])
			eps_fl = np.random.normal(0, proposal_stds[step_idx])
			proposed_depth_params[joint_depth_idx, idx_fp] += eps_fp
			proposed_depth_params[joint_depth_idx, idx_fl] += eps_fl

			if (proposed_depth_params[joint_depth_idx, idx_fp] < param_mins_depth[idx_fp, joint_depth_idx] or
				proposed_depth_params[joint_depth_idx, idx_fp] > param_maxs_depth[idx_fp, joint_depth_idx] or
				proposed_depth_params[joint_depth_idx, idx_fl] < param_mins_depth[idx_fl, joint_depth_idx] or
				proposed_depth_params[joint_depth_idx, idx_fl] > param_maxs_depth[idx_fl, joint_depth_idx]):
				continue_bounds = False

			if continue_bounds == True:
				fp_check, fl_check = _unconstrained_to_fractions(
					proposed_depth_params[joint_depth_idx, idx_fp],
					proposed_depth_params[joint_depth_idx, idx_fl])
				if fp_check < f_pyx_real_floor[joint_depth_idx] or fl_check < f_lherz_real_floor[joint_depth_idx]:
					continue_bounds = False

		if continue_bounds == True:

			object_backup = copy.deepcopy(object)

			if rand_dim < n_scalar_total:
				_lab_temp_proposed = proposed_scalars[1] if invert_lab_temp else lab_temp
				T_, P_, LAB = _update_geotherm(proposed_scalars[0], _lab_temp_proposed)
				object.set_temperature(T_)
				object.set_pressure(P_)

				#if bulk xFe is changed distributing iron among defined minerals.
				if 'bulk_xfe' in param_names:
					object.mantle_xfe_distribute(method = 'array')
				else:
					if triangle_calculation == True:
						for iz in range(n_depths):
							fp, fl = _unconstrained_to_fractions(
								proposed_depth_params[iz, idx_fp],
								proposed_depth_params[iz, idx_fl])
							object.f_pyx[iz] = fp
							object.f_lherz[iz] = fl

				if melt_thermodyn and melt_thermodyn_interp is not None:

					current_melt = np.zeros(n_depths)
					for iz in range(n_depths):
						T_C = object.T[iz] - 273.15
						water_wt = object.bulk_water[iz] * 1e-4
						melt_frac = float(melt_thermodyn_interp([T_C, water_wt, object.p[iz]]))
						if melt_frac < melt_frac_limit:
							melt_frac = 0.0
						current_melt[iz] = melt_frac
					object.melt_fluid_mass_frac[:n_depths] = current_melt

				# --- Set up fluid density interpolation for melt calculations ---
				if water_solv == True:

					object.mantle_water_distribute(method = 'array')
			elif rand_dim < n_scalar_total + n_param_total:

				#if one of the parameters are a composition parameter.
				if frac_bool[param_idx] == True:

					if triangle_calculation == True:
						raise AttributeError('Single compositional value and triangle calculation method cannot be used together.')

					if object.solid_phase_method == 2:
						_comp_list = [object.quartz_frac[depth_idx], object.plag_frac[depth_idx], object.amp_frac[depth_idx], object.kfelds_frac[depth_idx], object.opx_frac[depth_idx], object.cpx_frac[depth_idx],
							object.mica_frac[depth_idx], object.garnet_frac[depth_idx], object.sulphide_frac[depth_idx], object.graphite_frac[depth_idx], object.ol_frac[depth_idx], object.sp_frac[depth_idx], object.rwd_wds_frac[depth_idx],
							object.perov_frac[depth_idx], object.mixture_frac[depth_idx], object.other_frac[depth_idx]]
					else:
						_comp_list = [object.granite_frac[depth_idx],object.granulite_frac[depth_idx],object.sandstone_frac[depth_idx],object.gneiss_frac[depth_idx],object.amphibolite_frac[depth_idx],
							object.basalt_frac[depth_idx],object.mud_frac[depth_idx],object.gabbro_frac[depth_idx],object.other_rock_frac[depth_idx]]

					comp_old = _comp_list[comp_index[comp_index_sub]]

					comp_list = _comp_adjust_(np.array(_comp_list), proposed_depth_params[depth_idx, comp_index_sub], comp_old)

					for idx_t in range(len(_comp_list)):
						if object.solid_phase_method == 2:
							object.mineral_frac_list[idx_t][depth_idx] = comp_list[idx_t]
						else:
							object.rock_frac_list[idx_t][depth_idx] = comp_list[idx_t]

				#Changing the depth dependent parameter.
				if triangle_calculation == True and param_names[param_idx] in ('f_pyx', 'f_lherz'):
					# these columns hold logit-space values, so convert the PAIR
					# back to real fractions before handing them to the physics
					fp, fl = _unconstrained_to_fractions(
						proposed_depth_params[depth_idx, idx_fp],
						proposed_depth_params[depth_idx, idx_fl])
					object.f_pyx[depth_idx] = fp
					object.f_lherz[depth_idx] = fl
				else:
					getattr(object, param_names[param_idx])[depth_idx] = proposed_depth_params[depth_idx, param_idx]

				#if bulk xFe is changed distributing iron among defined minerals.
				if 'bulk_xfe' in param_names:
					object.mantle_xfe_distribute(method = 'index', sol_idx = depth_idx)

				if melt_thermodyn and melt_thermodyn_interp is not None:

					T_C = object.T[depth_idx] - 273.15
					water_wt = object.bulk_water[depth_idx] * 1e-4
					melt_frac = float(melt_thermodyn_interp([T_C, water_wt, object.p[depth_idx]]))
					if melt_frac < melt_frac_limit:
						melt_frac = 0.0
					object.melt_fluid_mass_frac[depth_idx] = melt_frac

				# --- Set up fluid density interpolation for melt calculations ---
				if water_solv == True:

					object.mantle_water_distribute(method = 'index', sol_idx = depth_idx)

			else:
				
				# joint f_pyx/f_lherz move: both columns were already updated
				# above (Site 4), just push the resulting real fractions onto
				# the object. Mirrors the single-move branch's post-processing
				# (melt/water recompute) for consistency, using its own
				# joint_depth_idx throughout -- never touches depth_idx.
				fp, fl = _unconstrained_to_fractions(
					proposed_depth_params[joint_depth_idx, idx_fp],
					proposed_depth_params[joint_depth_idx, idx_fl])
				object.f_pyx[joint_depth_idx] = fp
				object.f_lherz[joint_depth_idx] = fl

				if melt_thermodyn and melt_thermodyn_interp is not None:
					T_C = object.T[joint_depth_idx] - 273.15
					water_wt = object.bulk_water[joint_depth_idx] * 1e-4
					melt_frac = float(melt_thermodyn_interp([T_C, water_wt, object.p[joint_depth_idx]]))
					if melt_frac < melt_frac_limit:
						melt_frac = 0.0
					object.melt_fluid_mass_frac[joint_depth_idx] = melt_frac

				if water_solv == True:
					object.mantle_water_distribute(method = 'index', sol_idx = joint_depth_idx)

			#Calculating the conductivity and seismic velocities.
			if (vp_list is not None) or (vs_list is not None):
				v_bulk_, vp_, vs_ = object.calculate_seismic_velocities(method = 'array')

			#calculating conductivity later, because the conductivity may depend on gibbs-derived mineral
			#assemblage calculated within calculate_seismic_velocities call
			if cond_list is not None:
				cond_ = object.calculate_conductivity(method = 'array')

			if cond_list is not None:
				proposed_likelihood_cond, misf_cond = _likelihood(cond_, cond_list[index], sigma_cond[index])
				proposed_likelihood_cond = np.sum(proposed_likelihood_cond)
			else:
				proposed_likelihood_cond = 1
				misf_cond = np.zeros(len(object.T))

			if vp_list is not None:
				proposed_likelihood_vp, misf_vp = _likelihood(vp_, vp_list[index], sigma_vp[index], norm = 'linear')
				proposed_likelihood_vp = np.sum(proposed_likelihood_vp)
			else:
				proposed_likelihood_vp = 1
				misf_vp = np.zeros(len(object.T))
			if vs_list is not None:
				proposed_likelihood_vs, misf_vs = _likelihood(vs_, vs_list[index], sigma_vs[index], norm = 'linear')
				proposed_likelihood_vs = np.sum(proposed_likelihood_vs)
			else:
				proposed_likelihood_vs = 1
				misf_vs = np.zeros(len(object.T))
				
			if lab_depth is not None:
				proposed_likelihood_lab, misf_lab = _likelihood(LAB, lab_depth, sigma_lab, norm = 'linear')
			else:
				proposed_likelihood_lab = 1
				misf_lab = 0.0

			#Calculate prior likelihood for proposed parameters
			proposed_prior = 0.0
			if param_priors is not None:
				for ii in range(n_params):
					if param_priors[ii] is not None:
						prior_mean = param_priors[ii][0]   # array length n_depths
						prior_sigma = param_priors[ii][1]   # array length n_depths

						if triangle_calculation == True and ii == idx_fp:
							param_vals = np.array([_unconstrained_to_fractions(
								proposed_depth_params[iz, idx_fp], proposed_depth_params[iz, idx_fl])[0]
								for iz in range(n_depths)])
						elif triangle_calculation == True and ii == idx_fl:
							param_vals = np.array([_unconstrained_to_fractions(
								proposed_depth_params[iz, idx_fp], proposed_depth_params[iz, idx_fl])[1]
								for iz in range(n_depths)])
						else:
							param_vals = proposed_depth_params[:, ii]

						proposed_prior += np.sum(-0.5 * ((param_vals - prior_mean) / prior_sigma)**2)

			proposed_likelihood = np.exp(np.sum(misf_cond) + np.sum(misf_vp) + np.sum(misf_vs) + misf_lab + proposed_prior)

			if np.isnan(proposed_likelihood):
				n_reject_nan += 1
			else:
				# Calculate acceptance probability
				acceptance_ratio = proposed_likelihood / current_likelihood

				if np.random.rand() < acceptance_ratio:

					current_scalars = proposed_scalars.copy()
					current_depth_params = proposed_depth_params.copy()
					current_likelihood = proposed_likelihood

					if _ > burning:
						samples_SHF.append(current_scalars[0])
						if invert_lab_temp == True:
							samples_lab_temp.append(current_scalars[1])
						samples_lab.append(LAB)
						samples_temp.append(object.T.copy())
						samples_depth_params.append(_to_real_fractions(current_depth_params))
						misfits_cond.append(misf_cond)
						misfits_vp.append(misf_vp)
						misfits_vs.append(misf_vs)
						misfits_lab.append(misf_lab)
						if melt_thermodyn == True:
							melt_samples.append(object.melt_fluid_mass_frac[:n_depths].copy())
						for rname in record_names:
							samples_record[rname].append(np.array(getattr(object, rname)[:n_depths]).copy())
						accepted += 1
						n_accepted_per_dim[step_idx] += 1
				else:
					n_reject_likelihood += 1
					#reverting back to the old state.
					object = object_backup
					
		else:
			n_reject_bounds += 1

		if (_ - burning) > 0:
			acceptance_rate = accepted / (_ - burning)

		else:
			acceptance_rate = 0

		if _ > burning:
			# After base iterations, check if we need to continue
			if _ >= n_iter and (_ - n_iter) % 2000 == 0:
				if  ideal_acceptance_bounds[0] <= acceptance_rate <= ideal_acceptance_bounds[1]:
					print(f'Acceptance rate {acceptance_rate:.3f} is good. Terminating at {_} iterations.')
					break
				else:
					print(f'Acceptance rate {acceptance_rate:.3f} still not converged to desired acceptance rate. Continuing to maximum number of iterations...')

		acceptance_rates.append(acceptance_rate)
		misfits_all_cond.append(misf_cond.copy())
		misfits_all_vp.append(misf_vp.copy())
		misfits_all_vs.append(misf_vs.copy())
		misfits_all_lab.append(misf_lab.copy())
		samples_depth_params_all.append(_to_real_fractions(current_depth_params))
		samples_SHF_all.append(current_scalars[0])
		
		if invert_lab_temp == True:
			samples_lab_temp_all.append(current_scalars[1])
		samples_lab_all.append(LAB)
		samples_temp_all.append(object.T.copy())
		for rname in record_names:
			samples_record_all[rname].append(np.array(getattr(object, rname)[:n_depths]).copy())

		if melt_thermodyn == True:
			melt_samples_all.append(object.melt_fluid_mass_frac[:n_depths].copy())

		# Check if stuck after enough post-burn-in samples
		if _ != 0:
			if _ > burning and ((_ - burning) % 5000 == 0) and accepted == 0:
				print(text_color.RED + 'Zero acceptance after 5000 samples. Widening priors.' + text_color.END)
				if widen_count < max_widen_attempts:
					widen_count += 1
					# Widen priors for parameters that have them
					if param_priors is not None:
						for ii in range(n_params):
							if param_priors[ii] is not None:
								param_priors[ii][1][:] = param_priors[ii][1][:] * 1.25
								print(text_color.YELLOW + f'Index {index}: widening prior for {param_names[ii]} '
									f'to sigma_mean={param_priors[ii][1].mean():.3f}, attempt {widen_count}' + text_color.END)
		if acceptance_rate < 0.1:
			color_ovr = text_color.RED
		elif acceptance_rate < ideal_acceptance_bounds[0]:
			color_ovr = text_color.YELLOW
		elif acceptance_rate > ideal_acceptance_bounds[1]:
			color_ovr = text_color.YELLOW
		elif acceptance_rate > 0.5:
			color_ovr = text_color.RED
		else:
			color_ovr = text_color.GREEN

		if adaptive_alg == True:
			if (_ + 1) % adaptive_check_length == 0:
				print(f'Rejected based on likelihood per cent(%): {(1e2*n_reject_likelihood / _):.2f}')
				print(f'Rejected based on NaN petrophysical/thermodynamic calculation per cent(%): {(1e2*n_reject_nan / _):.2f}')
				print(f'Rejected based on bound limitations per cent(%): {(1e2*n_reject_bounds / _):.2f}')

				for d in range(n_step_dims):
					if n_attempted_per_dim[d] == 0:
						continue
					dim_acceptance = n_accepted_per_dim[d] / n_attempted_per_dim[d]
					if dim_acceptance < 0.1:
						proposal_stds[d] *= 0.8
						status, color = 'very low', text_color.RED
					elif dim_acceptance < ideal_acceptance_bounds[0]:
						proposal_stds[d] *= 0.95
						status, color = 'low', text_color.YELLOW
					elif dim_acceptance > 0.5:
						proposal_stds[d] *= 1.2
						status, color = 'very high', text_color.RED
					elif dim_acceptance > ideal_acceptance_bounds[1]:
						proposal_stds[d] *= 1.05
						status, color = 'high', text_color.YELLOW
					else:
						status, color = 'good', text_color.GREEN
					if step_size_limits is not None:
						proposal_stds[d] = min(proposal_stds[d], step_size_limits[d])

					if d < n_scalars:
						name = scalar_names[d]
					elif d < n_scalars + n_params:
						name = param_names[d - n_scalars]
					else:
						name = 'f_pyx_f_lherz_joint'
					print(color + f'  {name}: acceptance {status} ({dim_acceptance:.3f}), '
						f'attempted={n_attempted_per_dim[d]}, step={proposal_stds[d]:.4f}' + text_color.END)

				n_attempted_per_dim = [0] * n_step_dims
				n_accepted_per_dim = [0] * n_step_dims

				print(color_ovr + f'Overall acceptance: {acceptance_rate:.3f} | {(_/n_iter)*100:.1f}% done' + text_color.END)
		else:
			if (_ + 1) % adaptive_check_length == 0:
				print(text_color.GREEN + f'Acceptance Rate: {round(acceptance_rate,3)}' + text_color.END)


	misfits = [misfits_cond, misfits_vp, misfits_vs, misfits_lab]
	misfits_all = [misfits_all_cond, misfits_all_vp, misfits_all_vs, misfits_all_lab]

	samples_record = {k: np.array(v) for k, v in samples_record.items()}
	samples_record_all = {k: np.array(v) for k, v in samples_record_all.items()}

	results = {
		'samples_SHF': np.array(samples_SHF),
		'samples_lab': np.array(samples_lab),
		'samples_temp': np.array(samples_temp),
		'samples_depth_params': np.array(samples_depth_params),
		'acceptance_rates': np.array(acceptance_rates),
		'misfits': misfits,
		'samples_SHF_all': np.array(samples_SHF_all),
		'samples_lab_all': np.array(samples_lab_all),
		'samples_depth_params_all': np.array(samples_depth_params_all),
		'misfits_all': misfits_all,
		'samples_record': samples_record,
		'samples_record_all': samples_record_all,
	}
 
	if invert_lab_temp == True:
		results['samples_lab_temp'] = np.array(samples_lab_temp)
		results['samples_lab_temp_all'] = np.array(samples_lab_temp_all)
 
	if melt_thermodyn == True:
		results['melt_samples'] = np.array(melt_samples)
		results['melt_samples_all'] = np.array(melt_samples_all)
 
	return results
		
def _solv_MCMC_n_param(index, cond_list, object, initial_params, param_names, upper_limits,
	lower_limits, sigma_cond, proposal_stds, n_iter, burning, water_solv, comp_solv, melt_thermodyn, pres_interp, melt_frac_limit,
	vp_list = None, vs_list = None, sigma_vp = None, sigma_vs = None,
	adaptive_alg = True, ideal_acceptance_bounds = [0.2,0.3], adaptive_check_length = 1000,
	comp_index = [0,0], step_size_limits = None, transition_zone = False, param_priors = None,
	max_widen_attempts = 3, melt_thermodyn_interp = None):
	
	"""
	MCMC external solver for the metropolis_hastings_n_param function for parallelization purposes.
	Users should not call this function directly.
	"""
	
	widen_count = 0
	proposal_stds = list(proposal_stds)

	object = copy.deepcopy(object)
	
	if param_priors is not None:
		param_priors = copy.deepcopy(param_priors)

	if melt_thermodyn == True:
		
		if 'T' in param_names:
			idx_T = param_names.index('T')
		else:
			idx_T = None
		
		if 'bulk_water' in param_names:
			idx_B = param_names.index('bulk_water')
		else:
			idx_B = None

	n_params = len(param_names)
	#Using Metropolis-Hastings algorithm

	frac_bool = [False] * n_params
	comp_index_sub = None

	for ii in range(n_params):
		if 'frac' in param_names[ii]:
			if param_names[ii] != 'melt_fluid_mass_frac':
				frac_bool[ii] = True
				comp_index_sub = ii

	if sum(frac_bool) > 1:
		raise KeyError('Currently only one of the parameters chosen can be modal compositional parameter.')
		
	melt_solv = 'melt_fluid_mass_frac' in param_names

	current_params = np.array(initial_params[index])

	#Extracting per-index limits
	param_mins = np.array([lower_limits[i][index] for i in range(n_params)])
	param_maxs = np.array([upper_limits[i][index] for i in range(n_params)])

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

		comp_list = _comp_adjust_(np.array(_comp_list), current_params[comp_index_sub], comp_old)

		for idx_t in range(len(_comp_list)):
			if object.solid_phase_method == 2:
				object.mineral_frac_list[idx_t][index] = comp_list[idx_t]
			else:
				object.rock_frac_list[idx_t][index] = comp_list[idx_t]

		water_solv = True

	#Executing the commands - setting initial parameters
	for ii in range(n_params):
		getattr(object, param_names[ii])[index] = current_params[ii]
		
	if melt_thermodyn == True:

		if idx_T is not None:
			temp_ = current_params[idx_T] - 273.15
		else:
			temp_ = object.T[index] - 273.15
		if idx_B is not None:
			bw_ = current_params[idx_B] * 1e-4
		else:
			bw_ = object.bulk_water[index] * 1e-4

		if pres_interp == False:
			melt_frac = melt_thermodyn_interp([temp_,bw_])
		else:
			melt_frac = melt_thermodyn_interp([temp_,bw_, object.p[index]])

		if melt_frac < melt_frac_limit:
				
			melt_frac = np.array([0.0])

		getattr(object, 'melt_fluid_mass_frac')[index] = melt_frac

	if water_solv == True:

		if transition_zone == False:
			object.mantle_water_distribute(method = 'index', sol_idx = index)
		else:
			object.transition_zone_water_distribute(method = 'index', sol_idx = index)
			
		if (melt_solv == True) or (melt_thermodyn == True):
			
			#to interpolation of fluid density so eos do not have to be solved at each iteration.
			try:
				water_index = param_names.index('bulk_water')
				water_end = upper_limits[water_index]
				water_end = water_end[0]
			except ValueError:
				water_end = 1e5

			object.calculate_density_fluid(sol_idx = index, method = 'array', interp_for_iter = True, water_start = 0, water_end = water_end)

	if 'bulk_xfe' in param_names:

		object.mantle_xfe_distribute(method = 'index', sol_idx = index)

	#Calculating the initial conductivity
	cond_init = object.calculate_conductivity(method = 'index', sol_idx = index)
	if (vp_list is not None) or (vs_list is not None):
		v_bulk_init, vp_init, vs_init = object.calculate_seismic_velocities(method = 'index', sol_idx = index)

	current_likelihood_cond, current_misf = _likelihood(cond_init, cond_list[index], sigma_cond[index])
	
	if vp_list is not None:
		current_likelihood_vp, misf_vp = _likelihood(vp_init, vp_list[index], sigma_vp[index], norm = 'linear')
	else:
		current_likelihood_vp = 1
		misf_vp = 0
	if vs_list is not None:
		current_likelihood_vs, misf_vs = _likelihood(vs_init, vs_list[index], sigma_vs[index], norm = 'linear')
	else:
		current_likelihood_vs = 1
		misf_vs = 0

	#Calculate initial prior likelihood
	current_prior = 0.0
	if param_priors is not None:
		for ii in range(n_params):
			if param_priors[ii] is not None:
				prior_mean = param_priors[ii][0][index]
				prior_sigma = param_priors[ii][1][index]
				current_prior += -0.5 * ((current_params[ii] - prior_mean) / prior_sigma)**2

	current_likelihood = np.exp(current_misf + misf_vp + misf_vs + current_prior)

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
	melt_samples = []
	melt_samples_all = []
	accepted = 0
	print(text_color.GREEN + 'Monte-Carlo loop is started' + text_color.END)
	print(text_color.YELLOW + f'{n_iter*2} total samples, {n_iter} minimum samples.' + text_color.END)
	print(text_color.RED + f'{burning} burning samples.' + text_color.END)

	#loop for monte-carlo
	for _ in range(n_iter*2):
		
		#proposing the new parameters
		proposal = np.array(current_params)
		rand_node = np.random.randint(n_params) #random node generation 0 to n_params-1
		randomgen = np.random.normal(0, proposal_stds[rand_node], size=1)[0] #proposal for randomwalk step
		proposal[rand_node] = current_params[rand_node] + randomgen
		
		continue_bounds = np.all((proposal > param_mins) & (proposal < param_maxs))
		
		if continue_bounds == False:
			proposal = current_params #if out of bounds go back to the previous parameters

		if continue_bounds == True:
			#adjusting for composition if needed...

			object_backup = copy.deepcopy(object)

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

					comp_list = _comp_adjust_(np.array(_comp_list), proposal[comp_index_sub], comp_old)

					for idx_t in range(len(_comp_list)):
						if object.solid_phase_method == 2:
							object.mineral_frac_list[idx_t][index] = comp_list[idx_t]
						else:
							object.rock_frac_list[idx_t][index] = comp_list[idx_t]

			#setting up the proposed parameters
			for ii in range(n_params):
				getattr(object, param_names[ii])[index] = proposal[ii]
			
			if melt_thermodyn == True:

				if idx_T is not None:
					temp_ = proposal[idx_T] - 273.15
				else:
					temp_ = object.T[index] - 273.15
				if idx_B is not None:
					bw_ = proposal[idx_B] * 1e-4
				else:
					bw_ = object.bulk_water[index] * 1e-4

				if pres_interp == False:
					melt_frac = melt_thermodyn_interp([temp_,bw_])
				else:
					melt_frac = melt_thermodyn_interp([temp_,bw_, object.p[index]])
				
				if melt_frac < melt_frac_limit:
				
					melt_frac = np.array([0.0])

				getattr(object, 'melt_fluid_mass_frac')[index] = melt_frac

			#distribute water if needed
			if water_solv == True:
				if transition_zone == False:
					object.mantle_water_distribute(method = 'index', sol_idx = index)
				else:
					object.transition_zone_water_distribute(method = 'index', sol_idx = index)
			
			if 'bulk_xfe' in param_names:
		
				object.mantle_xfe_distribute(method = 'index', sol_idx = index)

			proposed_cond = object.calculate_conductivity(method = 'index', sol_idx = index)
			if (vp_list is not None) or (vs_list is not None):
				v_bulk, proposed_vp, proposed_vs = object.calculate_seismic_velocities(method = 'index', sol_idx = index)
			
			proposed_likelihood_cond, misf_cond = _likelihood(proposed_cond, cond_list[index], sigma_cond[index])
			
			if vp_list is not None:
				proposed_likelihood_vp, misf_vp = _likelihood(proposed_vp, vp_list[index], sigma_vp[index], norm = 'linear')
			else:
				proposed_likelihood_vp = 1
				misf_vp = 0
				
			if vs_list is not None:
				proposed_likelihood_vs, misf_vs = _likelihood(proposed_vs, vs_list[index], sigma_vs[index], norm = 'linear')
			else:
				proposed_likelihood_vs = 1
				misf_vs = 0
			
			#Calculate prior likelihood for proposed parameters
			proposed_prior = 0.0
			if param_priors is not None:
				for ii in range(n_params):
					if param_priors[ii] is not None:
						prior_mean = param_priors[ii][0][index]
						prior_sigma = param_priors[ii][1][index]
						proposed_prior += -0.5 * ((proposal[ii] - prior_mean) / prior_sigma)**2

			proposed_likelihood = np.exp(misf_cond + misf_vp + misf_vs + proposed_prior)

			# Calculate acceptance probability
			acceptance_ratio = proposed_likelihood / current_likelihood
			
			if np.random.rand() < acceptance_ratio:

				current_params = proposal
				current_likelihood = proposed_likelihood

				if _ > burning:
					samples.append(current_params.copy())
					misfits_cond.append(misf_cond)
					misfits_vp.append(misf_vp)
					misfits_vs.append(misf_vs)
					if melt_thermodyn == True:
						melt_samples.append(melt_frac.copy())
					accepted += 1

				
				if (_ - burning) > 0:
					acceptance_rate = accepted / (_ - burning)
					
				else:
					acceptance_rate = 0
					
				if _ > burning:
					# After base iterations, check if we need to continue
					if _ >= n_iter and (_ - n_iter) % 2000 == 0:
						if 0.2 <= acceptance_rate <= 0.3:
							print(f'Acceptance rate {acceptance_rate:.3f} is good. Terminating at {_} iterations.')
							break
						else:
							print(f'Acceptance rate {acceptance_rate:.3f} not converged to acceptable range. Continuing to maximum number of iterations...')
				acceptance_rates.append(acceptance_rate)
				misfits_all_cond.append(misf_cond)
				misfits_all_vp.append(misf_vp)
				misfits_all_vs.append(misf_vs)
				samples_all.append(current_params.copy())
				if melt_thermodyn == True:
					melt_samples_all.append(melt_frac.copy())
				
				# Check if stuck after enough post-burn-in samples
				if ((_ - burning) % 5000 == 0) and accepted == 0:
					print(text_color.RED + 'Zero acceptance after 5000 samples. Widening priors.' + text_color.END)
					if widen_count < max_widen_attempts:
						widen_count += 1
						# Widen priors for parameters that have them
						if param_priors is not None:
							for ii in range(n_params):
								if param_priors[ii] is not None:
									param_priors[ii][1][index] = param_priors[ii][1][index] * 1.25
									print(text_color.YELLOW + f'Index {index}: widening prior for {param_names[ii]} to sigma={param_priors[ii][1][index]:.1f}, attempt {widen_count}' + text_color.END)
				
				if adaptive_alg == True:
					if (_ + 1) % adaptive_check_length == 0:
						if acceptance_rate < 0.1:
							proposal_stds[rand_node] *= 0.8
							status = 'very low'
							color = text_color.RED
						elif acceptance_rate < ideal_acceptance_bounds[0]:
							proposal_stds[rand_node] *= 0.95
							status = 'low'
							color = text_color.YELLOW
						elif acceptance_rate > 0.5:
							proposal_stds[rand_node] *= 1.2
							status = 'very high'
							color = text_color.RED
						elif acceptance_rate > ideal_acceptance_bounds[1]:
							proposal_stds[rand_node] *= 1.05
							status = 'high'
							color = text_color.YELLOW
						else:
							status = 'good'
							color = text_color.GREEN
					
						# Enforce step size limits (both min and max)
						if step_size_limits is not None:
							if proposal_stds[rand_node] > step_size_limits[rand_node]:
								proposal_stds[rand_node] = step_size_limits[rand_node]
						print(color + f'Acceptance {status}: {acceptance_rate:.3f} | Steps: {[round(s,4) for s in proposal_stds]} | {(_/n_iter)*100:.1f}% done' + text_color.END)
				else:
					if (_ + 1) % adaptive_check_length == 0:
						print(text_color.GREEN + f'Acceptance Rate: {round(acceptance_rate,3)}' + text_color.END)
			else:
				#reverting back to the original object upon rejection of the sample.
				object = object_backup

	misfits = [misfits_cond, misfits_vp, misfits_vs]
	misfits_all = [misfits_all_cond, misfits_all_vp, misfits_all_vs]
	
	if melt_thermodyn == False:
		return np.array(samples), np.array(acceptance_rates), misfits, np.array(samples_all), np.array(misfits_all)
	else:
		return np.array(samples), np.array(acceptance_rates), misfits, np.array(samples_all), np.array(misfits_all), np.array(melt_samples), np.array(melt_samples_all)

def metropolis_hastings_n_param(object, cond_list, initial_params, param_names, upper_limits,
	lower_limits, sigma_cond, proposal_stds, n_iter, vp_list = None, vs_list = None, sigma_vs = None, sigma_vp = None,
	burning = 0, transition_zone = False, num_cpu = 1, param_priors = None, **kwargs):
 
	"""
	Perform Metropolis-Hastings MCMC inversion for electrical conductivity using n model parameters.
 
	This function uses a stochastic sampling approach to estimate the posterior distribution of n 
	input parameters based on observed conductivity and optionally seismic velocity data.
 
	Parameters
	----------
	object : object
		A pide model instance for calculating conductivity.
	cond_list : array-like
		Observed conductivity values to fit [S/m].
	initial_params : list or array-like
		Initial values for the n parameters to invert. Shape (n_points, n_params).
	param_names : list of str
		Names of the parameters to invert (must be attributes of `object`).
	upper_limits : tuple of array-like
		Upper bounds for each parameter. Shape (n_params,) where each element is array of length n_points.
	lower_limits : tuple of array-like
		Lower bounds for each parameter. Shape (n_params,) where each element is array of length n_points.
	sigma_cond : float or array-like
		Standard deviation or error for conductivity observations in logarithm of conductivity [S/m].
	proposal_stds : list or array-like
		Standard deviations for proposal distribution of the parameters. Length n_params.
	n_iter : int
		Number of iterations for the MCMC chain.
	vp_list : array-like, optional
		Observed Vp values to fit [km/s].
	vs_list : array-like, optional
		Observed Vs values to fit [km/s].
	sigma_vp : array-like, optional
		Standard deviation for Vp observations [km/s].
	sigma_vs : array-like, optional
		Standard deviation for Vs observations [km/s].
	burning : int, optional
		Number of initial iterations to discard as burn-in (default is 0).
	transition_zone : bool, optional
		If True, use transition zone water distribution functions (default is False).
	num_cpu : int, optional
		Number of CPU cores to use for parallel computation (default is 1).
	param_priors : list, optional
		List of length n_params. Each element is either None (flat prior) or a tuple 
		(prior_mean_array, prior_sigma_array) for a Gaussian prior. Arrays are length n_points.
		Example for 3 params (water=flat, melt=flat, T=Gaussian):
		param_priors = [None, None, (T_mean_array, sigma_T_array)]
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
		Bounds for the adaptive proposal step sizes. Length n_params.
	melt_thermodyn : bool, optional
		Method for determining melt from thermodynamic equations
	melt_frac_limit : float, optional
		Minimum melt fraction can be estimated by the inversion algorithm. Any value smaller than this value would be set to 0.
 
	Returns
	-------
	sample_distr : list of arrays
		Accepted parameter samples for each point.
	acceptance_rates : list of arrays
		Acceptance rate record for each point.
	misfits : list
		Misfits of the accepted distribution [cond, vp, vs].
	samples_all : list of arrays
		All parameter samples for each point.
	misfits_all : list
		All misfits [cond, vp, vs].
 
	Examples
	--------
	samples, acceptance_rates, misfits, samples_all, misfits_all = metropolis_hastings_n_param(
		object = p_obj, cond_list = [0.1, 0.1],
		initial_params = [[200, 0.25, 1300.0]],
		param_names = ['bulk_water', 'melt_fluid_mass_frac', 'T'],
		upper_limits = (np.array([2000,2000]), np.array([0.5,0.5]), np.array([1500,1500])),
		lower_limits = (np.array([0,0]), np.array([0,0]), np.array([1100,1100])),
		sigma_cond = [0.1, 0.1], proposal_stds = [200, 0.25, 50],
		n_iter = 2e5, burning = 1e4, transition_zone = False, num_cpu = 1,
		param_priors = [None, None, (T_mean_array, sigma_T_array)],
		adaptive_alg = True, step_size_limits = [25000, 0.5, 100])
	"""
 
	n_params = len(param_names)
 
	#Pre-checks for if
	if object.solid_phase_method == 2:
		object.set_mineral_water(ol = 0, opx = 0, cpx = 0, garnet = 0, mica = 0, amp = 0,
		quartz = 0, plag = 0, kfelds = 0, sulphide = 0, graphite = 0, sp = 0, rwd_wds = 0,
		perov = 0, mixture = 0, other = 0)
	elif object.solid_phase_method == 1:
		object.set_rock_water(granite = 0, granulite = 0, sandstone = 0, gneiss = 0,
		amphibolite = 0, basalt = 0, mud = 0, gabbro = 0, other_rock = 0)
 
	cond_check = object.calculate_conductivity()
 
	save_distr = kwargs.pop('save_distr', False)
	distr_file_names = kwargs.pop('distr_file_names', 'distribution_solution')
	adaptive_alg = kwargs.pop('adaptive_alg', True)
	ideal_acceptance_bounds = kwargs.pop('ideal_acceptance_bounds', [0.2, 0.3])
	adaptive_check_length = kwargs.pop('adaptive_check_length', 1000)
	step_size_limits = kwargs.pop('step_size_limits', None)
	melt_thermodyn = kwargs.pop('melt_thermodyn', False)
	melt_interp_object = kwargs.pop('melt_interp_object', None)
	melt_frac_limit = kwargs.pop('melt_frac_limit', 0.001)
 
	#Pre checks for the input parameters.
	if type(ideal_acceptance_bounds) == list:
		if len(ideal_acceptance_bounds) == 2:
			pass
		else:
			raise ValueError(f'ideal_acceptance_bounds has to be a list containing two values. Currently it is {ideal_acceptance_bounds}')
	else:
		raise ValueError(f'ideal_acceptance_bounds has to be a list containing two values. Currently it is {ideal_acceptance_bounds}')
 
	for name in param_names:
		try:
			getattr(object, name)
		except AttributeError:
			raise AttributeError(f'There is no such parameter name {name} for the pide object.')
  
	if burning >= n_iter:
		raise ValueError('Burning samples cannot be larger than the total iteration number (n_iter).')
 
	if len(cond_list) == len(initial_params) == len(upper_limits[0]) == len(lower_limits[0]) == len(sigma_cond) == len(object.T):
		pass
	else:
		raise IndexError('The length of the arrays for each conductivity solution (cond_list) are not same. cond_list, initial_params, upper_limits, lower_limits and sigma_conds has to be the same length.')
 
	#Check that all limits have the right number of parameters
	if len(upper_limits) != n_params:
		raise IndexError(f'upper_limits has {len(upper_limits)} entries but {n_params} parameters were specified.')
	if len(lower_limits) != n_params:
		raise IndexError(f'lower_limits has {len(lower_limits)} entries but {n_params} parameters were specified.')
	if len(proposal_stds) != n_params:
		raise IndexError(f'proposal_stds has {len(proposal_stds)} entries but {n_params} parameters were specified.')
	if len(initial_params[0]) != n_params:
		raise IndexError(f'initial_params has {len(initial_params[0])} entries but {n_params} parameters were specified.')
 
	min_list = ['quartz_frac', 'plag_frac', 'amp_frac', 'kfelds_frac', 'opx_frac', 'cpx_frac',
		'mica_frac', 'garnet_frac', 'sulphide_frac', 'graphite_frac', 'ol_frac', 'sp_frac', 'rwd_wds_frac',
		'perov_frac', 'mixture_frac', 'other_frac']
 
	rock_list = ['granite_frac', 'granulite_frac', 'sandstone_frac', 'gneiss_frac', 'amphibolite_frac',
			'basalt_frac', 'mud_frac', 'gabbro_frac', 'other_rock_frac']
 
	index_list = np.array(list(range(0, len(object.T)))) #creating the index array tied to the T array.
 
	comp_solv = False
	water_solv = False
	comp_type = None
	comp_index = []
	comp_type_list = []
 
	if any('water' in xx for xx in param_names) == True:
		if 'bulk_water' in param_names:
			water_solv = True
			if len(getattr(object, 'bulk_water')) != len(object.T):
				object.set_bulk_water(0.0)
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')
	
	
	if any('melt' in xx for xx in param_names) == True:
		if melt_thermodyn == False:
			water_solv = True
			for ii in range(n_params):
				if len(getattr(object, param_names[ii])) != len(object.T):
					object.set_parameter(param_names[ii], 0.0)
		else:
			raise KeyError('While melt_thermodyn is set to True, the user cannot choose melt fraction as a independent parameter. Melt fraction is estimated from thermodynamic equations.')
		
	if melt_thermodyn == True:
		water_solv = True
		if melt_interp_object is None:
			print('Establishing grid interpolator for thermodynamic melt modelling...')
			object.set_bulk_water(0.0)
			object.mantle_water_distribute()
			if np.all(object.p == object.p[0]):
				pres_interp = False
			else:
				pres_interp = True

			d_per_melt = (object.ol_frac * object.d_melt_ol) +\
					(object.opx_frac_wt * object.d_melt_opx) +\
					(object.cpx_frac_wt * object.d_melt_cpx) +\
					(object.garnet_frac_wt * object.d_melt_garnet)
			
			d_per_melt_avg = np.average(d_per_melt)

			if 'T' in param_names:
				idx_T = param_names.index('T')
				upper_lim_T = np.amax(upper_limits[idx_T]) - 273.15
				low_lim_T = np.amin(lower_limits[idx_T]) - 273.15
			else:
				upper_lim_T = np.amax(object.T) + 100.0 - 273.15
				low_lim_T = np.amin(object.T) - 273.15

			low_lim_X = 0
			up_lim_X = 5.0
			
			from pide.geodyn.mantlemelting.katz_2003 import F_wet
			
			T_grid = np.arange(low_lim_T, upper_lim_T, 25)
			
			if pres_interp == False:
				X_grid = np.linspace(low_lim_X, up_lim_X, 200)  # wt% water
				F_table = np.zeros((len(T_grid), len(X_grid)))
				for i, t in enumerate(T_grid):
					for j, x in enumerate(X_grid):
						F_table[i, j] = F_wet(T=t, P=object.p[0], X=x, D=d_per_melt_avg)
				F_table[F_table < 0] = 0.0

				melt_interp = RegularGridInterpolator((T_grid, X_grid), F_table, 
												bounds_error=False, fill_value=0.0)
				
			elif pres_interp == True:
				print('Startin thermodynamic pre-interpolator for melt. If this is taking more than 5 minutes restart.')
				low_lim_P = np.amin(object.p) - 0.1
				upper_lim_P = np.amax(object.p) + 0.1
				X_grid = np.linspace(low_lim_X, up_lim_X, 50)  # wt% water
				P_grid = np.linspace(low_lim_P,upper_lim_P,50)
				F_table = np.zeros((len(T_grid), len(X_grid), len(P_grid)))
				for i, t in enumerate(T_grid):
					for j, x in enumerate(X_grid):
						for k, pr in enumerate(P_grid):
							F_table[i, j, k] = F_wet(T=t, P=pr, X=x, D=d_per_melt_avg)

				F_table[F_table < 0] = 0.0
				print('F table finished')
				print('Interpolating')
				melt_interp = RegularGridInterpolator((T_grid, X_grid, P_grid), F_table, 
												bounds_error=False, fill_value=0.0)
				print('Interpolation table finished.')
		else:

			melt_interp = melt_interp_object
	else:

		melt_interp = None
		pres_interp = False
 
	if any('frac' in xx for xx in param_names) == True:

		comp_solv = True
		water_solv = True
 
		for ii in range(n_params):
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
				if len(getattr(object, param_names[ii])) != len(object.T):
					object.set_parameter(param_names[ii], 0.0)
 
		if (('mineral' in comp_type_list) == True) and (('rock' in comp_type_list) == True):
			raise ValueError('The user cannot enter both rock and mineral as the inversion parameter. Choose only one.')
 
	#The last check for setting up other parameters.
	for ii in range(n_params):
		if len(getattr(object, param_names[ii])) != len(object.T):
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
 
			process_item_partial = partial(_solv_MCMC_n_param, object = object, cond_list = cond_list,
			initial_params = initial_params, param_names = param_names,
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond,
			proposal_stds = proposal_stds, n_iter = n_iter, burning = burning,
			vp_list = vp_list, vs_list = vs_list, sigma_vp = sigma_vp, sigma_vs = sigma_vs,
			water_solv = water_solv, comp_solv = comp_solv, comp_index = comp_index,
			adaptive_alg = adaptive_alg,
			adaptive_check_length = adaptive_check_length, step_size_limits = step_size_limits,
			ideal_acceptance_bounds = ideal_acceptance_bounds, param_priors = param_priors,
			melt_thermodyn = melt_thermodyn, melt_thermodyn_interp = melt_interp, pres_interp = pres_interp,
			melt_frac_limit = melt_frac_limit)
 
			c = pool.map(process_item_partial, index_list)
 
		sample_distr = [x[0] for x in c]
		acceptance_rates = [x[1] for x in c]
		misfits = [x[2] for x in c]
		samples_all = [x[3] for x in c]
		misfits_all = [x[4] for x in c]
		if melt_thermodyn == True:
			melt_samples = [x[5] for x in c]
			melt_samples_all = [x[6] for x in c]
			
 
	else:
 
		sample_distr = []
		acceptance_rates = []
		misfits = []
		samples_all = []
		misfits_all = []
		melt_samples = []
		melt_samples_all = []
 
		for idx in range(0, len(index_list)):
			
			c = _solv_MCMC_n_param(index = index_list[idx], object = object, cond_list = cond_list,
			initial_params = initial_params, param_names = param_names,
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond,
			proposal_stds = proposal_stds, n_iter = n_iter, burning = burning,
			vp_list = vp_list, vs_list = vs_list, sigma_vp = sigma_vp, sigma_vs = sigma_vs,
			water_solv = water_solv, comp_solv = comp_solv, comp_index = comp_index,
			adaptive_alg = adaptive_alg,
			adaptive_check_length = adaptive_check_length, step_size_limits = step_size_limits,
			ideal_acceptance_bounds = ideal_acceptance_bounds, param_priors = param_priors,
			melt_thermodyn = melt_thermodyn, melt_thermodyn_interp = melt_interp, pres_interp = pres_interp,
			melt_frac_limit = melt_frac_limit)

			sample_distr.append(c[0])
			acceptance_rates.append(c[1])
			misfits.append(c[2])
			samples_all.append(c[3])
			misfits_all.append(c[4])
			if melt_thermodyn == True:
				melt_samples.append(c[5])
				melt_samples_all.append(c[6])
 
	if save_distr == True:
 
		from pide.utils.utils import save_h5_files
 
		array_names_idx = list(range(len(sample_distr)))
		array_names_idx = [str(element) for element in array_names_idx]
 
		save_h5_files(array_list=sample_distr, array_names=array_names_idx, file_name=distr_file_names + '_distr.h5')
		save_h5_files(array_list=acceptance_rates, array_names=array_names_idx, file_name=distr_file_names + '_acceptance.h5')
		save_h5_files(array_list=misfits, array_names=array_names_idx, file_name=distr_file_names + '_misfit.h5')
		save_h5_files(array_list=samples_all, array_names=array_names_idx, file_name=distr_file_names + '_distr_all.h5')
		save_h5_files(array_list=misfits_all, array_names=array_names_idx, file_name=distr_file_names + '_misfit_all.h5')
	
	if melt_thermodyn == False:
		return sample_distr, acceptance_rates, misfits, samples_all, misfits_all
	else:
		return sample_distr, acceptance_rates, misfits, samples_all, misfits_all, melt_samples, melt_samples_all
		
		
def metropolis_hastings_n_param_geotherm(object, depths, moho_depths,
	cond_list, vp_list, vs_list,
	sigma_cond_list, sigma_vp_list, sigma_vs_list,
	initial_SHF, initial_lab_temp,
	initial_params, param_names,
	upper_limits, lower_limits,
	SHF_bounds,
	proposal_stds, n_iter, burning, lab_depth = None, sigma_lab = None,
	melt_thermodyn=False,
	adaptive_alg=True, ideal_acceptance_bounds=[0.25, 0.45],
	adaptive_check_length=1000, step_size_limits=None,
	param_priors=None, SHF_prior=None, lab_temp = 1350,
	invert_lab_temp = False, lab_temp_bounds = None,
	composition = None, max_widen_attempts = 3,
	record_names = None,
	**kwargs):
 
	"""
	Perform Metropolis-Hastings MCMC inversion for electrical conductivity using n model parameters.
 
	This function uses a stochastic sampling approach to estimate the posterior distribution of n 
	input parameters based on observed conductivity and optionally seismic velocity data.
 
	Parameters
	----------
	object : object
		A pide model instance for calculating conductivity.
	cond_list : array-like
		Observed conductivity values to fit [S/m].
	initial_params : list or array-like
		Initial values for the n parameters to invert. Shape (n_points, n_params).
	param_names : list of str
		Names of the parameters to invert (must be attributes of `object`).
	upper_limits : tuple of array-like
		Upper bounds for each parameter. Shape (n_params,) where each element is array of length n_points.
	lower_limits : tuple of array-like
		Lower bounds for each parameter. Shape (n_params,) where each element is array of length n_points.
	sigma_cond : float or array-like
		Standard deviation or error for conductivity observations in logarithm of conductivity [S/m].
	proposal_stds : list or array-like
		Standard deviations for proposal distribution of the parameters. Length n_params.
	n_iter : int
		Number of iterations for the MCMC chain.
	vp_list : array-like, optional
		Observed Vp values to fit [km/s].
	vs_list : array-like, optional
		Observed Vs values to fit [km/s].
	sigma_vp : array-like, optional
		Standard deviation for Vp observations [km/s].
	sigma_vs : array-like, optional
		Standard deviation for Vs observations [km/s].
	burning : int, optional
		Number of initial iterations to discard as burn-in (default is 0).
	transition_zone : bool, optional
		If True, use transition zone water distribution functions (default is False).
	num_cpu : int, optional
		Number of CPU cores to use for parallel computation (default is 1).
	param_priors : list, optional
		List of length n_params. Each element is either None (flat prior) or a tuple 
		(prior_mean_array, prior_sigma_array) for a Gaussian prior. Arrays are length n_points.
		Example for 3 params (water=flat, melt=flat, T=Gaussian):
		param_priors = [None, None, (T_mean_array, sigma_T_array)]
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
		Bounds for the adaptive proposal step sizes. Length n_params.
	melt_thermodyn : bool, optional
		Method for determining melt from thermodynamic equations
	melt_frac_limit : float, optional
		Minimum melt fraction can be estimated by the inversion algorithm. Any value smaller than this value would be set to 0.
 
	Returns
	-------
	sample_distr : list of arrays
		Accepted parameter samples for each point.
	acceptance_rates : list of arrays
		Acceptance rate record for each point.
	misfits : list
		Misfits of the accepted distribution [cond, vp, vs].
	samples_all : list of arrays
		All parameter samples for each point.
	misfits_all : list
		All misfits [cond, vp, vs].
 
	Examples
	--------
	samples, acceptance_rates, misfits, samples_all, misfits_all = metropolis_hastings_n_param(
		object = p_obj, cond_list = [0.1, 0.1],
		initial_params = [[200, 0.25, 1300.0]],
		param_names = ['bulk_water', 'melt_fluid_mass_frac', 'T'],
		upper_limits = (np.array([2000,2000]), np.array([0.5,0.5]), np.array([1500,1500])),
		lower_limits = (np.array([0,0]), np.array([0,0]), np.array([1100,1100])),
		sigma_cond = [0.1, 0.1], proposal_stds = [200, 0.25, 50],
		n_iter = 2e5, burning = 1e4, transition_zone = False, num_cpu = 1,
		param_priors = [None, None, (T_mean_array, sigma_T_array)],
		adaptive_alg = True, step_size_limits = [25000, 0.5, 100])
	"""
 
	n_params = len(param_names)
 
	#Pre-checks for if
	if object.solid_phase_method == 2:
		object.set_mineral_water(ol = 0, opx = 0, cpx = 0, garnet = 0, mica = 0, amp = 0,
		quartz = 0, plag = 0, kfelds = 0, sulphide = 0, graphite = 0, sp = 0, rwd_wds = 0,
		perov = 0, mixture = 0, other = 0)
	elif object.solid_phase_method == 1:
		object.set_rock_water(granite = 0, granulite = 0, sandstone = 0, gneiss = 0,
		amphibolite = 0, basalt = 0, mud = 0, gabbro = 0, other_rock = 0)
 
	cond_check = object.calculate_conductivity()
 
	save_distr = kwargs.pop('save_distr', False)
	distr_file_names = kwargs.pop('distr_file_names', 'distribution_solution')
	adaptive_alg = kwargs.pop('adaptive_alg', True)
	ideal_acceptance_bounds = kwargs.pop('ideal_acceptance_bounds', [0.2, 0.3])
	adaptive_check_length = kwargs.pop('adaptive_check_length', 1000)
	step_size_limits = kwargs.pop('step_size_limits', None)
	melt_thermodyn = kwargs.pop('melt_thermodyn', False)
	melt_interp_object = kwargs.pop('melt_interp_object', None)
	melt_frac_limit = kwargs.pop('melt_frac_limit', 0.001)
 
	#Pre checks for the input parameters.
	if type(ideal_acceptance_bounds) == list:
		if len(ideal_acceptance_bounds) == 2:
			pass
		else:
			raise ValueError(f'ideal_acceptance_bounds has to be a list containing two values. Currently it is {ideal_acceptance_bounds}')
	else:
		raise ValueError(f'ideal_acceptance_bounds has to be a list containing two values. Currently it is {ideal_acceptance_bounds}')
 
	for name in param_names:
		try:
			getattr(object, name)
		except AttributeError:
			raise AttributeError(f'There is no such parameter name {name} for the pide object.')
  
	if burning >= n_iter:
		raise ValueError('Burning samples cannot be larger than the total iteration number (n_iter).')
 
	if len(cond_list) == len(initial_params) == len(upper_limits[0]) == len(lower_limits[0]) == len(sigma_cond) == len(object.T):
		pass
	else:
		raise IndexError('The length of the arrays for each conductivity solution (cond_list) are not same. cond_list, initial_params, upper_limits, lower_limits and sigma_conds has to be the same length.')
 
	#Check that all limits have the right number of parameters
	if len(upper_limits) != n_params:
		raise IndexError(f'upper_limits has {len(upper_limits)} entries but {n_params} parameters were specified.')
	if len(lower_limits) != n_params:
		raise IndexError(f'lower_limits has {len(lower_limits)} entries but {n_params} parameters were specified.')
	if len(proposal_stds) != n_params:
		raise IndexError(f'proposal_stds has {len(proposal_stds)} entries but {n_params} parameters were specified.')
	if len(initial_params[0]) != n_params:
		raise IndexError(f'initial_params has {len(initial_params[0])} entries but {n_params} parameters were specified.')
 
	min_list = ['quartz_frac', 'plag_frac', 'amp_frac', 'kfelds_frac', 'opx_frac', 'cpx_frac',
		'mica_frac', 'garnet_frac', 'sulphide_frac', 'graphite_frac', 'ol_frac', 'sp_frac', 'rwd_wds_frac',
		'perov_frac', 'mixture_frac', 'other_frac']
 
	rock_list = ['granite_frac', 'granulite_frac', 'sandstone_frac', 'gneiss_frac', 'amphibolite_frac',
			'basalt_frac', 'mud_frac', 'gabbro_frac', 'other_rock_frac']
 
	index_list = np.array(list(range(0, len(object.T)))) #creating the index array tied to the T array.
 
	comp_solv = False
	water_solv = False
	comp_type = None
	comp_index = []
	comp_type_list = []
 
	if any('water' in xx for xx in param_names) == True:
		if 'bulk_water' in param_names:
			water_solv = True
			if len(getattr(object, 'bulk_water')) != len(object.T):
				object.set_bulk_water(0.0)
		else:
			raise ValueError('You cannot change just a single phase water content. If you are after fitting for a single phase, try bulk_water as the parameter.')
	
	
	if any('melt' in xx for xx in param_names) == True:
		if melt_thermodyn == False:
			water_solv = True
			for ii in range(n_params):
				if len(getattr(object, param_names[ii])) != len(object.T):
					object.set_parameter(param_names[ii], 0.0)
		else:
			raise KeyError('While melt_thermodyn is set to True, the user cannot choose melt fraction as a independent parameter. Melt fraction is estimated from thermodynamic equations.')
		
	if melt_thermodyn == True:
		water_solv = True
		if melt_interp_object is None:
			print('Establishing grid interpolator for thermodynamic melt modelling...')
			object.set_bulk_water(0.0)
			object.mantle_water_distribute()
			if np.all(object.p == object.p[0]):
				pres_interp = False
			else:
				pres_interp = True

			d_per_melt = (object.ol_frac * object.d_melt_ol) +\
					(object.opx_frac_wt * object.d_melt_opx) +\
					(object.cpx_frac_wt * object.d_melt_cpx) +\
					(object.garnet_frac_wt * object.d_melt_garnet)
			
			d_per_melt_avg = np.average(d_per_melt)

			if 'T' in param_names:
				idx_T = param_names.index('T')
				upper_lim_T = np.amax(upper_limits[idx_T]) - 273.15
				low_lim_T = np.amin(lower_limits[idx_T]) - 273.15
			else:
				upper_lim_T = np.amax(object.T) + 100.0 - 273.15
				low_lim_T = np.amin(object.T) - 273.15

			low_lim_X = 0
			up_lim_X = 5.0
			
			from pide.geodyn.mantlemelting.katz_2003 import F_wet
			
			T_grid = np.arange(low_lim_T, upper_lim_T, 25)
			
			if pres_interp == False:
				X_grid = np.linspace(low_lim_X, up_lim_X, 200)  # wt% water
				F_table = np.zeros((len(T_grid), len(X_grid)))
				for i, t in enumerate(T_grid):
					for j, x in enumerate(X_grid):
						F_table[i, j] = F_wet(T=t, P=object.p[0], X=x, D=d_per_melt_avg)
				F_table[F_table < 0] = 0.0

				melt_interp = RegularGridInterpolator((T_grid, X_grid), F_table, 
												bounds_error=False, fill_value=0.0)
				
			elif pres_interp == True:
				print('Startin thermodynamic pre-interpolator for melt. If this is taking more than 5 minutes restart.')
				low_lim_P = np.amin(object.p) - 0.1
				upper_lim_P = np.amax(object.p) + 0.1
				X_grid = np.linspace(low_lim_X, up_lim_X, 50)  # wt% water
				P_grid = np.linspace(low_lim_P,upper_lim_P,50)
				F_table = np.zeros((len(T_grid), len(X_grid), len(P_grid)))
				for i, t in enumerate(T_grid):
					for j, x in enumerate(X_grid):
						for k, pr in enumerate(P_grid):
							F_table[i, j, k] = F_wet(T=t, P=pr, X=x, D=d_per_melt_avg)

				F_table[F_table < 0] = 0.0
				print('F table finished')
				print('Interpolating')
				melt_interp = RegularGridInterpolator((T_grid, X_grid, P_grid), F_table, 
												bounds_error=False, fill_value=0.0)
				print('Interpolation table finished.')
		else:

			melt_interp = melt_interp_object
	else:

		melt_interp = None
		pres_interp = False

	if ('f_pyx' in param_names) or ('f_lherz' in param_names):
		triangle_calculation = True
	else:
		triangle_calculation = False
 
	if any('frac' in xx for xx in param_names) == True:

		comp_solv = True
		water_solv = True
 
		for ii in range(n_params):
			if param_names[ii] != 'melt_fluid_mass_frac':
				if triangle_calculation == True:
					raise AttributeError('f_pyx or f_lherz can not be considered at the same time as modal /'
					'mineralogy fractions as inversion paramteres..')
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
				if len(getattr(object, param_names[ii])) != len(object.T):
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
 
			process_item_partial = partial(_solv_MCMC_column, object = object, cond_list = cond_list,
			initial_params = initial_params, param_names = param_names,
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond,
			proposal_stds = proposal_stds, n_iter = n_iter, burning = burning,
			vp_list = vp_list, vs_list = vs_list, sigma_vp = sigma_vp, sigma_vs = sigma_vs,
			comp_index = comp_index,
			adaptive_alg = adaptive_alg,
			adaptive_check_length = adaptive_check_length, step_size_limits = step_size_limits,
			ideal_acceptance_bounds = ideal_acceptance_bounds, param_priors = param_priors,
			melt_thermodyn = melt_thermodyn, melt_thermodyn_interp = melt_interp, pres_interp = pres_interp,
			melt_frac_limit = melt_frac_limit)
 
			c_ = pool.map(process_item_partial, index_list)

			
 
	else:
 
		sample_distr = []
		acceptance_rates = []
		misfits = []
		samples_all = []
		misfits_all = []
		melt_samples = []
		melt_samples_all = []
 
		for idx in range(0, len(index_list)):
			
			c_ = _solv_MCMC_n_param(index = index_list[idx], object = object, cond_list = cond_list,
			initial_params = initial_params, param_names = param_names,
			upper_limits = upper_limits, lower_limits = lower_limits, sigma_cond = sigma_cond,
			proposal_stds = proposal_stds, n_iter = n_iter, burning = burning,
			vp_list = vp_list, vs_list = vs_list, sigma_vp = sigma_vp, sigma_vs = sigma_vs,
			comp_index = comp_index,
			adaptive_alg = adaptive_alg,
			adaptive_check_length = adaptive_check_length, step_size_limits = step_size_limits,
			ideal_acceptance_bounds = ideal_acceptance_bounds, param_priors = param_priors,
			melt_thermodyn = melt_thermodyn, melt_thermodyn_interp = melt_interp, pres_interp = pres_interp,
			melt_frac_limit = melt_frac_limit)

			sample_distr.append(c_[0])
			acceptance_rates.append(c_[1])
			misfits.append(c_[2])
			samples_all.append(c_[3])
			misfits_all.append(c_[4])
			if melt_thermodyn == True:
				melt_samples.append(c_[5])
				melt_samples_all.append(c_[6])
 
	

	array_names_idx = list(range(len(sample_distr)))
	array_names_idx = [str(element) for element in array_names_idx]

	save_h5_files(array_list=sample_distr, array_names=array_names_idx, file_name=distr_file_names + '_distr.h5')
	save_h5_files(array_list=acceptance_rates, array_names=array_names_idx, file_name=distr_file_names + '_acceptance.h5')
	save_h5_files(array_list=misfits, array_names=array_names_idx, file_name=distr_file_names + '_misfit.h5')
	save_h5_files(array_list=samples_all, array_names=array_names_idx, file_name=distr_file_names + '_distr_all.h5')
	save_h5_files(array_list=misfits_all, array_names=array_names_idx, file_name=distr_file_names + '_misfit_all.h5')
	