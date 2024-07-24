#!/usr/bin/env python3

import numpy as np
import sys

def interpolate_2d_fields(mesh_field, vals, mesh_out, method = 'linear'):

	"""
	To interpolate 2d fields
	
	mesh_field: input mesh in np.meshgrid
	vals: Value field associated with the mesh_field
	mesh_out: output mesh to interpolate data to
	method: method of interpolations see scipy.interpolate.griddata to choose between options
	
	"""

	#Function to interpolate a larger data array defined by np.meshgrid(mesh_field) and associated values(vals)

	from scipy.interpolate import griddata
	
	
	points_x = []
	points_y = []
	for i in range(0,len(mesh_field[0])):
		for j in range(0,len(mesh_field[0][0])):	
			points_x.append(mesh_field[0][i][j])
			points_y.append(mesh_field[1][i][j])
			
	points_interp = np.column_stack((np.array(points_x), np.array(points_y)))
		
	interp_vals = griddata(points_interp, vals, (mesh_out[0],mesh_out[1]),method = method)
	
	return interp_vals
		
def interpolate_3d_fields(mesh_tuple, vals, mesh_out, method = 'linear'):

	"""
	To interpolate 3d fields using scipy.RegularGridInterpolator
	
	mesh_field: input mesh in tuple (xmesh,y_mesh,z_mesh)
	vals: 1D Value field associated with the mesh_field
	mesh_out: output mesh to interpolate data to in tupl (xmesh,y_mesh,z_mesh), mostly to use for move to mesh_centers
	
	"""
	
	from scipy.interpolate import RegularGridInterpolator as rgi
	
	x_len = len(mesh_tuple[0])
	y_len = len(mesh_tuple[1])
	z_len = len(mesh_tuple[2])
	
	#checking if z is ascending, to overcome scipy 'strictly ascending' requirement.
	#the outputs should be the same since we only deal with dimensions in this space.
	if mesh_tuple[2][1] < mesh_tuple[2][0]:
		ascending_bool = True
	else:
		ascending_bool = False
		
	if ascending_bool == True:
	
		mesh_tuple_inv = (mesh_tuple[0],mesh_tuple[1], -1 * mesh_tuple[2])
		mesh_out_inv = (mesh_out[0],mesh_out[1], -1 * mesh_out[2])
		
	else:
	
		mesh_tuple_inv = mesh_tuple
		mesh_out_inv = mesh_out
	
	matrix = np.zeros((x_len,y_len,z_len))
	try:
		#filling the matrix with 1-d array values based on the xyz dimensions
		a = 0

		for k in range(0,z_len):
			for j in range(0,y_len):
				for i in range(0,x_len):
					if a != len(vals):
						matrix[i][j][k] = vals[a]
					a = a + 1
					
	except IndexError as e:
		
		import traceback
		
		traceback.print_exc()
		print('pideErrorHelp: There is a mismatch between the entered value array and mesh parameters. This likely results from the mesh indices are not associated with the value array used.')
		sys.exit()
	
	interp_func = rgi(mesh_tuple_inv, matrix, method = method) #interpolation function in 3-D space.
		
	interp_array = []
	for k in range(0,len(mesh_out_inv[2])):
		for j in range(0,len(mesh_out_inv[1])):
			for i in range(0,len(mesh_out_inv[0])):
		
				interp_array.append(np.array([mesh_out_inv[0][i],mesh_out_inv[1][j],mesh_out_inv[2][k]]))
			
	interp_array = np.array(interp_array)
		
	interpolated_vals = interp_func(interp_array) #getting out the interpolated values
	
	return interpolated_vals
	