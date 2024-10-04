#!/usr/bin/env python3

import numpy as np

def return_material_bool(material_index,model_array, material_skip, model_type):
	
	"""A function to get indexes of a pide.Material object of a material array given in 
	pide.Model object. 	
	"""
	
	array_bool = np.where(model_array == material_index)
	
	if model_type == "underworld_2d":
		
		if material_skip is not None:
			# Recreate arrays with every material_skip from each original array
			_new_array_bool = tuple(arr[::material_skip] for arr in array_bool)
			
		else:
		
			_new_array_bool = array_bool
			
	elif model_type == "underworld_3d":
		
		array_bool = array_bool[0]#tuple is not neccesary
		if material_skip != None:
			# Recreate arrays with every material_skip from each original array
			_new_array_bool = array_bool[::material_skip]
			_new_array_bool = np.array(_new_array_bool)
			
		else:
			_new_array_bool = array_bool
			
	else:
	
		raise ValueError('The model_type has entered wrongly.')
	
	return _new_array_bool
		

