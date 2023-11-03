#!/usr/bin/env python3

import os, csv
import numpy as np

def convert_2DModel_2_MARE2DEM(file_out, conductivity_array, mesh, boundaries = None, **kwargs):

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
		
	#converting mesh in kilometers into meters.
	if mesh_unit == 'kilometres':
		mesh[0] = mesh[0] * 1e3
		mesh[1] = mesh[1] * 1e3
	elif mesh_unit == 'metres':
		pass
	else:
		raise ValueError('Please enter a valid response for mesh_unit: "kilometres" or "metres".')
		
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
				
	#finding mesh_center
	mc_lateral = (boundaries['right'] - boundaries['left']) / 2.0
	
	#centering the mesh and boundaries with the center point in [y] direction.
	mesh[0] = mesh[0] - mc_lateral
	boundaries['right'] = boundaries['right'] - mc_lateral
	boundaries['left'] = boundaries['left'] - mc_lateral
	
	#appending lines to write from te conductivity array and mesh
	lines = []

	for i in range(0,len(conductivity_array)):
		for j in range(0,len(conductivity_array[0])):
			
			if (mesh[0][i][j] <= boundaries['right']) and (mesh[0][i][j] >= boundaries['left']):
				if (mesh[1][i][j] <= boundaries['bottom']) and (mesh[1][i][j] >= boundaries['top']):
					lines.append(','.join((str(mesh[0][i][j]),str(mesh[1][i][j]),str(conductivity_array[i][j][0]) + '\n')))
	
	#saving the file
	filesave_composition = open(file_out ,'w')
	filesave_composition.writelines(lines)
	filesave_composition.close()
	
	print('The synthetic conductivity model had been succesfully written as MARE2DEM/Mamba input at: ')
	print(file_out)