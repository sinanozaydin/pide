#!/usr/bin/env python3

import SEL
import numpy as np

def return_material_bool(material_index,model_array, material_skip):

	#getting the material indexes
	array_bool = np.where(model_array == material_index)
	
	
	if material_skip != None:
		# Recreate arrays with every material_skip from each original array
		_new_array_bool = tuple(arr[::material_skip] for arr in array_bool)
		
	else:
	
		_new_array_bool = array_bool
	
	return _new_array_bool
	

