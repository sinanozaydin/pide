#!/usr/bin/env python3

import sys, csv
import numpy as np
import h5py

def read_h5_files(temp_h5, pressure_h5, mesh_h5, material_h5, strain_h5, stress_h5, melt_h5, strain_rate_h5):

	try:
		strain_file = h5py.File(strain_h5, 'r')
	except FileNotFoundError:
		print('The strain rate file is not found.')

	try:
		temp_file = h5py.File(temp_h5, 'r')

	except FileNotFoundError:
		print('The temperature file is not found.')
		
	try:
		pressure_file = h5py.File(pressure_h5, 'r')

	except FileNotFoundError:
		print('The pressure file is not found.')
		
	try:
		mesh_file = h5py.File(mesh_h5, 'r')

	except FileNotFoundError:
		print('The mesh file is not found.')

	try:
		material_file = h5py.File(material_h5, 'r')

	except FileNotFoundError:
		print('The material file is not found.')
		
	try:
		stress_file = h5py.File(stress_h5, 'r')

	except FileNotFoundError:
		print('The stress file is not found.')
		
	try:
		melt_file = h5py.File(melt_h5, 'r')

	except FileNotFoundError:
		print('The melt file is not found.')
		
	try:
		strain_rate_file = h5py.File(strain_rate_h5, 'r')

	except FileNotFoundError:
		print('The strain rate file is not found.')
		
	
	return temp_file, pressure_file, mesh_file, material_file, strain_file, stress_file, melt_file, strain_rate_file
	
def read_uw_material_names_from_py_input(py_file):

	#Simple function for reading csv files and give out filtered output for given delimiter (delim)

	file_obj = open(py_file,'rt',encoding = "utf8") #Creating file object
	file_csv = csv.reader(file_obj) #Reading the file object with csv module, delimiter assigned to ','
	py_file_text = [] #Creating empty array to append data

	#Appending data from csb object
	for row in file_csv:
		py_file_text.append(row)

	#Filtering data for None elements read.
	for j in range(0,len(py_file_text)):
		py_file_text[j] = list(filter(None,py_file_text[j]))
	py_file_text = list(filter(None,py_file_text))

	
	idx_material_py = []

	for i in range(0,len(py_file_text)):
		if ('.add_material' in py_file_text[i][0]) == True:
			idx_material_py.append(i)
	
	material_names = []

	for i in idx_material_py:
		idx_local = py_file_text[i][0].find('=')
		material_names.append(py_file_text[i][0][:idx_local].strip())
	
	return material_names
	
def setup_2d_mesh(mesh_data):

	mesh_array = np.array(list(mesh_data['vertices']))
	mesh_x_array = mesh_array[:,0]
	mesh_y_array = -1 * mesh_array[:,1] #changing the minus direction in the Earth.

	print('#######################')
	print('Setting up the mesh parameters...')
	for i in range(0,3):
		print('                    ')
	#setting up mesh parameters
	
	#finding max and miny

	max_x = np.amax(mesh_x_array)
	min_x = np.amin(mesh_x_array)

	max_y = np.amax(mesh_y_array)
	min_y = np.amin(mesh_y_array)
	
	borders_mesh = [max_x, min_x, max_y, min_y]

	print('Maximum X:  ' + str(max_x) + '   km')
	print('Minimum X:  ' + str(min_x) + '   km')
	print('Maximum Y:  ' + str(min_y) + '   km')
	print('Minimum Y:  ' + str(min_y) + '   km')

	increment_in_x = np.abs(mesh_x_array[1] - mesh_x_array[0])
	x_steps = int((max_x - min_x) / increment_in_x)
	
	increment_in_y = np.abs(mesh_y_array[x_steps+1] - mesh_y_array[0])
	y_steps = int((max_y - min_y) / increment_in_y)

	x_mesh = np.arange(min_x, max_x + increment_in_x, increment_in_x)
	y_mesh = np.arange(max_y, min_y - increment_in_y, -increment_in_y)

	x_mesh_centers = x_mesh[:-1] + (increment_in_x / 2.0)
	y_mesh_centers = y_mesh[:-1] - (increment_in_y / 2.0)

	mesh = np.meshgrid(x_mesh, y_mesh)
	mesh_center = np.meshgrid(x_mesh_centers, y_mesh_centers)
	
	return mesh, mesh_center, x_mesh, y_mesh, borders_mesh
		
def setup_material(material_data, material_names):

	material_array = np.array(list(material_data['data']))

	for i in range(0,len(material_array)):
		material_array[i] = round(material_array[i][0]) #rounding to get rid of that areas with 

	print(' ')
	print('Materials included in the py startup file, matching up with the projMaterial.h5 material index identifiers.')
	unique_materials = np.unique(material_array)
	air_material_idx = []
	print('id    materialname')
	for i in range(0,len(unique_materials)):

		print(str(int(unique_materials[i])) + '    ' + str(material_names[i]))

		if (('air' in material_names[i]) or ('Air' in material_names[i])) == True:
			air_material_idx.append(i)
			
	return material_array, air_material_idx
	
def setup_uw_data_array_PROJ(data):

	array = np.array(list(data['data']))
	
	return array
	
def plot_2D_underworld_Field(x_array = None, y_array = None, Field = None,cblimit_up = None, cblimit_down = None, log_bool = False,cb_name = 'coolwarm'):

	import matplotlib.pyplot as plt
	import matplotlib.colors as colors

	fig = plt.figure(figsize = (12,7))
	ax = plt.subplot(111)
	ax.set_ylim(np.amax(y_array),np.amin(y_array))
	ax.set_xlim(np.amin(x_array),np.amax(x_array))
	if log_bool == True:
		cax = ax.scatter(x_array, y_array ,c = Field, cmap = cb_name, norm=colors.LogNorm(), marker = 's', linewidth = 0.005, edgecolor = 'k')
	elif log_bool == False:
		cax = ax.scatter(x_array, y_array ,c = Field, cmap = cb_name, marker = 's', linewidth = 0.005, edgecolor = 'k')
	cax.set_clim(cblimit_down,cblimit_up)

	if log_bool == True:
		bondary = np.logspace(np.log10(cblimit_down),np.log10(cblimit_up))
		tick_array = np.arange(np.log10(cblimit_down),np.log10(cblimit_up)+1, 1)
		tick_array_list = 10.0**tick_array
		cbar_cax = fig.colorbar(cax,boundaries=bondary ,orientation="vertical", pad=0.05,
		ticks = tick_array_list, ax = ax)
		
	elif log_bool == False:
		bondary = np.linspace(cblimit_down, cblimit_up)
		cbar_cax = fig.colorbar(cax,boundaries=bondary ,orientation="vertical", pad=0.05, ax = ax)
		
		
	plt.show()