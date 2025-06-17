import numpy as np
import csv
from scipy.optimize import minimize

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
	
def read_csv(filename,delim,linefiltering = True):
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
		raise ValueError("Sum of new fixed values must be less than 100")
	
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
	return results
	
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
