import numpy as np

def check_type(input):

	if isinstance(input,(int,float)):
		method_calc = 'scalar'
	else:
		method_calc = 'array'
		
	return method_calc