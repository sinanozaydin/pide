import numpy as np

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