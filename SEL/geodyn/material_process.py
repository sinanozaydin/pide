#!/usr/bin/env python3

import SEL
import numpy as np

def return_material_bool(material_index,model_array):

	array_bool = np.where(model_array == material_index)
			
	return array_bool
	

