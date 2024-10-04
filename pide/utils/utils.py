import numpy as np
import csv

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
	
def read_csv(filename,delim):
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
			
			process_item_partial = partial(_associate_coordinates_, x_target = target_x,
			y_target = target_y, x_sample = sample_x, y_sample = sample_y)
			
			c = pool.map(process_item_partial, index_list)
			
		idx_array = [x for x in c]
			
	else:
	
		idx_array = []
		
		for i in index_list:
		
			idx = _associate_coordinates_(i, x_target = target_x,
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
