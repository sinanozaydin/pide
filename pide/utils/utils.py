import numpy as np
import csv
from scipy.optimize import minimize
import xarray as xr

def _associate_coordinates_(index, x_target, y_target, x_sample, y_sample):
	
	"""
	Function to associate_coordinates method for parallelisation purposes.
	
	Input:
	int: index - index for the x_target array
	array: x_target - array for the target search x-direction
	array: y_target -array for the target search in y-direction
	array: x_sample -array for the samples in x-direction
	array: y_sample -array for the samples in y-direction
	
	Output:
	int: idx_final - the closest index where x_sample and y_sample is closes to x_target[index] y_target[index]
	"""

	idx_target = (np.abs(y_sample-y_target[index])).argmin()
	idx_target_lists = [idx for idx, value in enumerate(y_sample) if value == y_sample[idx_target]]
	idx_ = (np.abs(x_sample[idx_target_lists]-x_target[index])).argmin()
	idx_final = idx_target_lists[idx_]
	
	return idx_final

def _associated_coordinates_2(index, x_target, y_target, x_sample, y_sample):

	distance = np.sqrt((y_target[index] - y_sample)**2.0 + ((x_target[index] - x_sample)**2.0))
	idx_final = np.argmin(distance)

	return idx_final

def check_type(input):

	"""
	A check type function
	"""

	if isinstance(input,(int,float,np.float64, np.int64)):
		tip = 'scalar'
	elif isinstance(input, str):
		tip = 'string'
	elif isinstance(input, dict):
		tip = 'dict'
	else:
		tip = 'array'
		
	return tip
	
def array_modifier(input, array, varname):

	"""
	A function to extend the input into the length of the given array. In pide this is mostly used for
	arrays to match with temperature array (self.T).
	"""
		
	if type(input) == int:
		
		ret_array = np.ones(len(array)) * input
		
	elif type(input) == float:
		
		ret_array = np.ones(len(array)) * input
		
	elif type(input) == np.float64:
		
		ret_array = np.ones(len(array)) * input
		
	elif type(input) == list:
		
		if len(input) == 1:
			ret_array = np.ones(len(array)) * input[0]
		else:
			ret_array = np.array(input)
			if len(ret_array) != len(array):
				
				raise RuntimeError('The entered list of ***' + varname + '*** does not match the length of the entered temperature array.')
		
	elif type(input) == np.ndarray:
	
		if len(input) == 1:
			ret_array = np.ones(len(array)) * input[0]
		else:
			ret_array = input
			if len(ret_array) != len(array):
				raise RuntimeError('The entered list of ***' + varname + '*** does not match the length of the entered temperature array.')
			
	return ret_array
	
def read_csv(filename,delim,linefiltering = True,elementfiltering = True):
	"""
	Simple function for reading csv files and give out filtered output for given delimiter (delim)
	
	Input:
	str: filename - filename string or full path to the csv file.
	str: delim - delimiter for the csv file.
	
	Output:
	array: data
	"""
	
	with open(filename,'rt',encoding = "utf8") as file_obj:
	
		file_csv = csv.reader(file_obj,delimiter = delim) #Reading the file object with csv module, delimiter assigned to ','
		data = [] #Creating empty array to append data
	
		#Appending data from csv object
		for row in file_csv:
			data.append(row)
	
		#Filtering data for None elements read.
		if elementfiltering == True:
			for j in range(0,len(data)):
				data[j] = list(filter(None,data[j]))
			
		if linefiltering == True:
			data = list(filter(None,data))
	
		return data
		
		
def associate_coordinates(sample_x, sample_y, target_x, target_y,  num_cpu = 1, filename = 'idx.mat' ,method = 'return'):

	"""
	Function to associate_coordinates method for parallelisation purposes.
	
	Input:
	array: x_target - array for the target search x-direction
	array: y_target -array for the target search in y-direction
	array: x_sample -array for the samples in x-direction
	array: y_sample -array for the samples in y-direction
	
	Output:
	array: idx_array - the closest index where x_sample and y_sample is closes to x_target y_target
	"""

	if num_cpu != 1:
	
		import multiprocessing
		from functools import partial
	
	index_list = range(0,len(target_x))
	
	if num_cpu != 1:
	
		with multiprocessing.Pool(processes=num_cpu) as pool:
			
			process_item_partial = partial(_associated_coordinates_2, x_target = target_x,
			y_target = target_y, x_sample = sample_x, y_sample = sample_y)
			
			c = pool.map(process_item_partial, index_list)
			
		idx_array = [x for x in c]
			
	else:
	
		idx_array = []
		
		for i in index_list:
		
			idx = _associated_coordinates_2(i, x_target = target_x,
			y_target = target_y, x_sample = sample_x, y_sample = sample_y) 
			idx_array.append(idx)
		
	
	idx_array = np.array(idx_array)
	
	if method == 'save':
	
		from scipy.io import savemat
		savemat(filename, {'matrix': idx_array})
		print('File is saved at the location: ' + str(filename))
		
	else:
	
		return idx_array
	
def sort_through_external_list(first_list, second_list):
	# A function to sort secondary list dependent on the sorting of the first list.
	combined_lists = zip(first_list, second_list)
	sorted_combined_lists = sorted(combined_lists)
	sorted_second_list = [element[1] for element in sorted_combined_lists]
	return sorted_second_list
	
def save_h5_files(array_list, array_names, file_name = "Data.h5"):

	import h5py
	
	with h5py.File(file_name, 'w') as f:
		
		for i in range(len(array_list)):
		
			f.create_dataset(array_names[i], data = array_list[i])

	print(f'Results are saved at: {file_name}')
	
def save_mcmc_results_nc(filename, lon, lat, sample_distr, acceptance_rates, misfits_all,
	param_names, depth=None, extra_data=None, melt_thermodyn=False, melt_samples=None,
	title='MCMC petrophysical inversion results'):

	n_points = len(lon)
	n_params = len(param_names)

	# Calculate statistics from sample distributions
	mean_params = np.full((n_points, n_params), np.nan)
	std_params = np.full((n_points, n_params), np.nan)
	median_params = np.full((n_points, n_params), np.nan)
	p05_params = np.full((n_points, n_params), np.nan)
	p95_params = np.full((n_points, n_params), np.nan)
	final_acceptance = np.full(n_points, np.nan)
	mean_misfit_cond = np.full(n_points, np.nan)
	mean_misfit_vp = np.full(n_points, np.nan)
	mean_misfit_vs = np.full(n_points, np.nan)

	# Melt statistics if thermodynamic melt is used
	if melt_thermodyn == True and melt_samples is not None:
		melt_mean = np.full(n_points, np.nan)
		melt_std = np.full(n_points, np.nan)
		melt_median = np.full(n_points, np.nan)
		melt_p05 = np.full(n_points, np.nan)
		melt_p95 = np.full(n_points, np.nan)

	for i in range(n_points):

		if sample_distr[i] is not None and len(sample_distr[i]) > 0:
			if sample_distr[i][0] is not None:
				samples = np.array(sample_distr[i])
				if samples.ndim == 2 and samples.shape[0] > 1:
					mean_params[i, :] = np.mean(samples, axis=0)
					std_params[i, :] = np.std(samples, axis=0)
					median_params[i, :] = np.median(samples, axis=0)
					p05_params[i, :] = np.percentile(samples, 5, axis=0)
					p95_params[i, :] = np.percentile(samples, 95, axis=0)

		if acceptance_rates[i] is not None and len(acceptance_rates[i]) > 0:
			if acceptance_rates[i][-1] is not None:
				final_acceptance[i] = acceptance_rates[i][-1]

		if misfits_all[i] is not None:
			try:
				if misfits_all[i][0] is not None and len(misfits_all[i][0]) > 0:
					mean_misfit_cond[i] = np.mean(misfits_all[i][0])
				if misfits_all[i][1] is not None and len(misfits_all[i][1]) > 0:
					mean_misfit_vp[i] = np.mean(misfits_all[i][1])
				if misfits_all[i][2] is not None and len(misfits_all[i][2]) > 0:
					mean_misfit_vs[i] = np.mean(misfits_all[i][2])
			except (TypeError, IndexError):
				pass

		if melt_thermodyn == True and melt_samples is not None:
			if melt_samples[i] is not None and len(melt_samples[i]) > 0:
				m = np.array(melt_samples[i]).flatten()
				if len(m) > 1:
					melt_mean[i] = np.mean(m)
					melt_std[i] = np.std(m)
					melt_median[i] = np.median(m)
					melt_p05[i] = np.percentile(m, 5)
					melt_p95[i] = np.percentile(m, 95)

	# Build the dataset
	data_vars = {}

	# Add parameter statistics
	for j in range(n_params):
		name = param_names[j]
		data_vars[f'{name}_mean'] = (['point'], mean_params[:, j])
		data_vars[f'{name}_std'] = (['point'], std_params[:, j])
		data_vars[f'{name}_median'] = (['point'], median_params[:, j])
		data_vars[f'{name}_p05'] = (['point'], p05_params[:, j])
		data_vars[f'{name}_p95'] = (['point'], p95_params[:, j])

	# Add thermodynamic melt statistics
	if melt_thermodyn == True and melt_samples is not None:
		data_vars['melt_thermodyn_mean'] = (['point'], melt_mean)
		data_vars['melt_thermodyn_std'] = (['point'], melt_std)
		data_vars['melt_thermodyn_median'] = (['point'], melt_median)
		data_vars['melt_thermodyn_p05'] = (['point'], melt_p05)
		data_vars['melt_thermodyn_p95'] = (['point'], melt_p95)

	# Add diagnostics
	data_vars['acceptance_rate'] = (['point'], final_acceptance)
	data_vars['mean_misfit_cond'] = (['point'], mean_misfit_cond)
	data_vars['mean_misfit_vp'] = (['point'], mean_misfit_vp)
	data_vars['mean_misfit_vs'] = (['point'], mean_misfit_vs)

	# Add extra data if provided
	if extra_data is not None:
		for key, value in extra_data.items():
			if isinstance(value, tuple):
				data_vars[key] = (['point'], value[0])
			else:
				data_vars[key] = (['point'], value)

	# Create dataset
	dataset = xr.Dataset(
		data_vars,
		coords={
			'lon': (['point'], np.array(lon)),
			'lat': (['point'], np.array(lat)),
		}
	)

	# Add depth as attribute if provided
	if depth is not None:
		dataset.attrs['depth_km'] = depth

	dataset.attrs['title'] = title
	dataset.attrs['param_names'] = str(param_names)
	dataset.attrs['melt_thermodyn'] = str(melt_thermodyn)

	# Add coordinate attributes
	dataset['lon'].attrs = {'units': 'degrees_east', 'long_name': 'Longitude'}
	dataset['lat'].attrs = {'units': 'degrees_north', 'long_name': 'Latitude'}

	# Add parameter attributes
	param_units = {
		'bulk_water': 'ppm',
		'melt_fluid_mass_frac': 'fraction',
		'T': 'C',
		'ol_frac': 'fraction',
		'opx_frac': 'fraction',
		'cpx_frac': 'fraction',
		'garnet_frac': 'fraction',
	}

	for j in range(n_params):
		name = param_names[j]
		units = param_units.get(name, '')
		for suffix in ['_mean', '_std', '_median', '_p05', '_p95']:
			if f'{name}{suffix}' in dataset:
				dataset[f'{name}{suffix}'].attrs = {'units': units, 'long_name': f'{name} {suffix[1:]}'}

	# Add melt thermodynamic attributes
	if melt_thermodyn == True and melt_samples is not None:
		for suffix in ['_mean', '_std', '_median', '_p05', '_p95']:
			dataset[f'melt_thermodyn{suffix}'].attrs = {
				'units': 'fraction',
				'long_name': f'Thermodynamic melt fraction {suffix[1:]}'
			}

	dataset['acceptance_rate'].attrs = {'units': '', 'long_name': 'Final MCMC acceptance rate'}
	dataset['mean_misfit_cond'].attrs = {'units': 'log10(S/m)', 'long_name': 'Mean conductivity misfit'}
	dataset['mean_misfit_vp'].attrs = {'units': 'km/s', 'long_name': 'Mean Vp misfit'}
	dataset['mean_misfit_vs'].attrs = {'units': 'km/s', 'long_name': 'Mean Vs misfit'}

	# Add extra data attributes
	if extra_data is not None:
		for key, value in extra_data.items():
			if isinstance(value, tuple) and len(value) >= 2:
				attrs = {'units': value[1]}
				if len(value) >= 3:
					attrs['long_name'] = value[2]
				dataset[key].attrs = attrs

	# Write to file
	dataset.to_netcdf(filename)
	print(f'Saved MCMC results to {filename}')
	print(f'  {n_points} points, {n_params} parameters: {param_names}')
	if melt_thermodyn:
		print(f'  Includes thermodynamic melt fraction statistics')
	print(f'  Variables: {list(dataset.data_vars)}')

	return dataset
	
def _comp_adjust_idx_based(_comp_list, comp_alien, idx, array = False):

	"""A method to adjust composition of one mineral/rock without considering the replacement weights
	"""
	
	if array == False:
		ratio = (comp_alien - _comp_list[idx]) / (np.sum(_comp_list) - _comp_list[idx])
		comp_list = _comp_list - (_comp_list * ratio)
		comp_list[idx] = comp_alien
	else:
		ratio = (comp_alien - _comp_list[:,idx]) / (np.sum(_comp_list, axis = 1) - _comp_list[:,idx])
		comp_list = _comp_list - (_comp_list.T * ratio).T
		comp_list[:,idx] = comp_alien
	
	return comp_list

def _all_equal(arrays):

	"""
	Return True if every sub-list in `arrays` is equal to the first one.
	"""
	if not arrays:
		return True   # empty → trivially “all the same”
	first = arrays[0]
	return all(sub == first for sub in arrays)


def _comp_adjust_melts(sio2,na2o,k2o,comp_dict_rest):
	
	# Assume wt% total must be 100
	# SiO2 and Na2O are given, rest are variables
	def objective(x):  # x = [Al2O3, FeO, MgO, CaO, K2O, TiO2]
		total = sum(x) + sio2 + na2o + k2o
		return abs(total - 100)

	initial_guess = [comp_dict_rest['Al2O3'], comp_dict_rest['MgO'], comp_dict_rest['FeO'], comp_dict_rest['CaO'],
					comp_dict_rest['TiO2'],comp_dict_rest['MnO'], comp_dict_rest['P2O5'], 
					comp_dict_rest['Cr2O3']]  # make educated guess
	bounds = [(0, 25)] * 6  # basic bounds

	result = minimize(objective, initial_guess, bounds=bounds)

	comp_adjusted = [sio2, result.x[0], result.x[1], result.x[2],
				  result.x[3], na2o, k2o, result.x[4],result.x[5],
				  result.x[6],result.x[7], 0.0]

	return comp_adjusted

def modify_melt_composition(composition, indexes_to_change, new_values):
	"""
	Adjust composition array(s) by changing specific values while maintaining sum = 100.
	
	Parameters:
	composition: Can be either:
		- Single composition: list/array of numbers
		- Multiple compositions: list/array of lists/arrays
	indexes_to_change: list of int - indexes of values to change (same for all compositions)
	new_values: Can be either:
		- For single composition: list of numbers
		- For multiple compositions: list of lists (one per composition)
	
	Returns:
	- For single composition: numpy array
	- For multiple compositions: list of numpy arrays
	"""
	# Detect if we have multiple compositions
	if isinstance(composition[0], (list, np.ndarray)):
		# Multiple compositions case
		return adjust_composition_batch(composition, indexes_to_change, new_values)
	else:
		# Single composition case
		return adjust_single_composition(composition, indexes_to_change, new_values)

def adjust_single_composition(composition, indexes_to_change, new_values):
	"""
	Adjust a single composition array by changing specific values while maintaining sum = 100.
	
	Parameters:
	composition: list or numpy array - the original composition (should sum to ~100)
	indexes_to_change: list of int - indexes of values to change
	new_values: list of float - new values for the specified indexes
	
	Returns:
	numpy array - adjusted composition that sums to 100
	"""
	# Convert to numpy array for easier manipulation
	comp = np.array(composition, dtype=float)
	
	# Validate inputs
	if len(indexes_to_change) != len(new_values):
		raise ValueError("Number of indexes must match number of new values")
	
	if any(idx >= len(comp) for idx in indexes_to_change):
		raise ValueError("Index out of range")
	
	# Calculate sum of new fixed values
	sum_new_fixed = sum(new_values)
	
	if sum_new_fixed >= 100:
		raise ValueError("Components of the melt composition makes up more than 100 percent. This usually occurs"+\
		"when a user accidentally sets the melt water content more than 1e5 ppm (100%wt.). Try to set a lower bound for bulk water content in your inversion.")
	
	# Create a mask for values that will NOT be changed
	mask = np.ones(len(comp), dtype=bool)
	mask[indexes_to_change] = False
	
	# Get current sum of values that will remain flexible
	current_flexible_sum = comp[mask].sum()
	
	# Calculate target sum for flexible values
	target_flexible_sum = 100 - sum_new_fixed
	
	# Handle edge case where no flexible values remain
	if current_flexible_sum == 0:
		if target_flexible_sum > 0:
			raise ValueError("Cannot redistribute to zero-sum flexible values")
	else:
		# Calculate scaling factor for flexible values
		scaling_factor = target_flexible_sum / current_flexible_sum
		
		# Apply scaling to flexible values
		comp[mask] *= scaling_factor
	
	# Set the new fixed values
	for idx, new_val in zip(indexes_to_change, new_values):
		comp[idx] = new_val
	
	return comp

def adjust_composition_batch(compositions, indexes_to_change, new_values_per_composition):
	"""
	Adjust multiple composition arrays.
	
	Parameters:
	compositions: list of lists/arrays - multiple compositions
	indexes_to_change: list of int - indexes to change (same for all compositions)
	new_values_per_composition: list of lists - new values for each composition
	
	Returns:
	list of numpy arrays - adjusted compositions
	"""
	if len(compositions) != len(new_values_per_composition):
		raise ValueError("Number of compositions must match number of new value arrays")
	
	results = []
	for i, comp in enumerate(compositions):
		adjusted = adjust_single_composition(comp, indexes_to_change, new_values_per_composition[i])
		results.append(adjusted)
	return np.array(results)
	
class text_color:
   
   #color object for to be called by the print outs.
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'
