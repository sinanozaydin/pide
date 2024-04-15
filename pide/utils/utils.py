import numpy as np
import csv

def _associate_coordinates_(index, x_target, y_target, x_sample, y_sample):

	idx_target = (np.abs(y_sample-y_target[index])).argmin()
	idx_target_lists = [idx for idx, value in enumerate(y_sample) if value == y_sample[idx_target]]
	idx_ = (np.abs(x_sample[idx_target_lists]-x_target[index])).argmin()
	idx_final = idx_target_lists[idx_]
	
	return idx_final

def check_type(input):

	if isinstance(input,(int,float)):
		method_calc = 'scalar'
	else:
		method_calc = 'array'
		
	return method_calc
	
def array_modifier(input, array, varname):
		
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

	#Simple function for reading csv files and give out filtered output for given delimiter (delim)
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
		
		
def associate_coordinates(sample_x, sample_y, target_x, target_y, filename, num_cpu = 1, method = 'return'):

	"""
	Returns the associated indexes 
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
	
def check_return(func):
	def wrapper(*args, **kwargs):
		result = func(*args, **kwargs)
		if not wrapper.var_assigned:
			return None
		return result

	wrapper.var_assigned = False
	return wrapper