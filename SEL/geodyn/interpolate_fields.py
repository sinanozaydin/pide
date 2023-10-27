#!/usr/bin/env python3

import numpy as np

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
