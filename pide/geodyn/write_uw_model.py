#!/usr/bin/env python3

import numpy as np
import h5py

def write_2d_field_h5(Field, filename_out, nan_placeholder = -999, nan_interpolate = False, xmesh = None, ymesh = None):

	if nan_interpolate == True:
				
		Field = np.squeeze(Field, axis=1)
		
		contains_nan = np.isnan(Field).any()
		
		if contains_nan == True:
		
			from scipy.interpolate import griddata
			#interpolating 
			
			xi = xmesh
			yi = ymesh
			x_i,  y_i = np.meshgrid(xmesh,ymesh)
			
			mask = np.isnan(Field)
			points = np.column_stack((x_i[~mask], y_i[~mask]))
			values = Field[~mask]
					
			Field = griddata(points, values, (x_i, y_i), method = 'linear')

	data_2_write = []

	for i in range(0,len(Field)):
	
		for j in range(0,len(Field[0])):
			
			if np.isnan(Field[i][j]) == True:
				data_2_write.append(nan_placeholder)
			else:
				data_2_write.append(Field[i][j])

	with h5py.File(filename_out, 'w') as file:
	
		file.create_dataset('data', data = data_2_write)
		
	print('The file has written as:' + str(filename_out))
	
def write_3d_field_h5(Field, filename_out, nan_placeholder = -999, nan_interpolate = False, xmesh = None, ymesh = None, zmmesh = None):

	data_2_write = []

	for i in range(0,len(Field)):
	
		if np.isnan(Field[i]) == True:
			data_2_write.append(nan_placeholder)
		else:
			data_2_write.append(Field[i])

	with h5py.File(filename_out, 'w') as file:
	
		file.create_dataset('data', data = data_2_write)
		
	print('The file has written as:' + str(filename_out))
		
		