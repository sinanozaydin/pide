#!/usr/bin/env python3

import os, csv
import numpy as np

def convert_2DModel_2_MARE2DEM(file_out, conductivity_array, mesh, boundaries = None, **kwargs):

	import scipy.io

	#This is a function that converts the constructed 2D geodynamic model into data points
	#that can be read by Mamba software that is used for creating model files for MARE2DEM.
	"""
	Input Paramaters
	-----------------
	file_out: filename or full directory to output the result. Includes the format (e.g., csv)
	conductivity_array: calculated conductivity field in np.ndarray()
	mesh: mesh of the associated conductivity field.
	boundaries: boundaries dictionary in  {'left','right','top','bottom'}	
	
	cond_unit: 'conductivity' or 'resistivity'
	mesh_unit: 'kilometres' or 'metres'
	"""
	
	#optional kwarg calls
	cond_unit = kwargs.pop('cond_unit', 'conductivity')
	mesh_unit = kwargs.pop('mesh_unit', 'kilometres')
	
	if cond_unit == 'conductivity':
		conductivity_array = 1.0 / conductivity_array
	elif cond_unit == 'resistivity':
		pass
	else:
		raise ValueError('Please enter a valid response for cond_unit: "conductivity" or "resistivity".')
		
	
		
	boundary_list = ['left','right','top','bottom']
	#Defaulting into mesh boundaries if w
	if boundaries == None:
		boundaries = {'left': np.amin(mesh[0]), 'right': np.amax(mesh[0]), 'top': 0, 'bottom':np.amax(mesh[1])}
	else:
		for item in boundary_list:
			if item not in boundaries:
				if item == 'top':
					boundaries[item] = 0.0
				elif item == 'right':
					boundaries[item] = np.amax(mesh[0])
				elif item == 'left':
					boundaries[item] = np.amin(mesh[0])
				elif item == 'bottom':
					boundaries[item] = np.amax(mesh[1])
						
	#converting mesh in kilometers into meters.
	if mesh_unit == 'kilometres':
		mesh[0] = mesh[0] * 1e3
		mesh[1] = mesh[1] * 1e3
		boundaries['top'] = boundaries['top'] * 1e3
		boundaries['right'] = boundaries['right'] * 1e3
		boundaries['left'] = boundaries['left'] * 1e3
		boundaries['bottom'] = boundaries['bottom'] * 1e3
	elif mesh_unit == 'metres':
		pass
	else:
		raise ValueError('Please enter a valid response for mesh_unit: "kilometres" or "metres".')
		
	#finding mesh_center
	mc_lateral = (boundaries['right'] - boundaries['left']) / 2.0
	
	#centering the mesh and boundaries with the center point in [y] direction.
	mesh[0] = mesh[0] - mc_lateral
	boundaries['right'] = boundaries['right'] - mc_lateral
	boundaries['left'] = boundaries['left'] - mc_lateral
	
	#appending lines to write from te conductivity array and mesh
	lines = []
	array_y = []
	array_z = []
	array_rho = []
	array_save = []
	
	for i in range(0,len(conductivity_array)):
		for j in range(0,len(conductivity_array[0])):
			if np.isnan(conductivity_array[i][j][0]) == False:
				if (mesh[0][i][j] <= boundaries['right']):
					if (mesh[0][i][j] >= boundaries['left']):
						if (mesh[1][i][j] <= boundaries['bottom']):
							if (mesh[1][i][j] >= boundaries['top']):
								# array_save.append([mesh[0][i][j],mesh[1][i][j],conductivity_array[i][j][0]])
								array_y.append(mesh[0][i][j])
								array_z.append(mesh[1][i][j])
								array_rho.append(conductivity_array[i][j][0])
								# lines.append('  '.join((str(mesh[0][i][j]),str(mesh[1][i][j]),str(np.log10(conductivity_array[i][j][0])) + '\n')))
	
	#saving the file
	# filesave_composition = open(file_out ,'w')
	# filesave_composition.writelines(lines)
	# filesave_composition.close()
	data_2_write = {
	'y': array_y,
	'z': array_z,
	'rho': array_rho}
	
	# data_2_write = {
	# 'yzrho': array_save
	# }
	
	scipy.io.savemat(file_out, data_2_write)
	
	print('The synthetic conductivity model had been succesfully written as MARE2DEM/Mamba input at: ')
	print(file_out)
	print('##################')
	print('Model bounds entered:')
	print('--------- ' + str(boundaries['top']) + '---------')
	print('|                                |')
	print('|                                |')
	print(str(boundaries['left']) + '                    ' + str(boundaries['right']))
	print('|                                |')
	print('|                                |')
	print('--------- ' + str(boundaries['bottom']) + '---------')