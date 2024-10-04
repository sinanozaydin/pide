#!/usr/bin/env python3

import csv
import numpy as np
import h5py

def read_h5_file(h5):

	try:
		h5_file = h5py.File(h5, 'r')
	except FileNotFoundError:
		print('The entered file is not found')

	return h5_file

def read_uw_material_names_from_py_input(py_file):

	"""A function to scrape material names and indexes from a underworld python input.
	Input:
	str: py_file - full path or filename to the python file.
	Output:
	Materil name lists.
	"""

	#Simple function for reading csv files and give out filtered output for given delimiter (delim)

	with open(py_file,'rt',encoding = "utf8") as file_obj:
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

	"""
	A function to convert 2D mesh_data to a usable array.

	Input Parameters:
	mesh_data read from the file

	Output parameters:
	mesh in np.meshgrid
	mesh_center in np.meshgrid
	x_mesh: mesh array in x-dir in np.array
	y_mesh: mesh array in y-dir in np.array
	x_mesh_centers: mesh centers array in x-dir in np.array
	y_mesh_centers: mesh centers array in y-dir in np.array
	"""

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
	print('Maximum Y:  ' + str(max_y) + '   km')
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

	return mesh, mesh_center, x_mesh, y_mesh, x_mesh_centers, y_mesh_centers, borders_mesh

def setup_3d_mesh(mesh_data):

	"""
	A function to convert 3D mesh_data to a usable array.

	Input Parameters:
	mesh_data read from the file

	Output parameters:
	mesh in np.meshgrid
	mesh_center in np.meshgrid
	x_mesh: mesh array in x-dir in np.array
	y_mesh: mesh array in y-dir in np.array
	z_mesh: mesh array in z-dir in np.array
	x_mesh_centers: mesh centers array in x-dir in np.array
	y_mesh_centers: mesh centers array in y-dir in np.array
	z_mesh_centers: mesh centers array in z-dir in np.array
	"""

	mesh_array = np.array(list(mesh_data['vertices']))
	mesh_x_array = mesh_array[:,0]
	mesh_y_array = mesh_array[:,1]
	mesh_z_array = -1 * mesh_array[:,2] #changing the minus direction in the Earth.
	#finding max and miny

	max_x = np.amax(mesh_x_array)
	min_x = np.amin(mesh_x_array)

	max_y = np.amax(mesh_y_array)
	min_y = np.amin(mesh_y_array)

	max_z = np.amax(mesh_z_array)
	min_z = np.amin(mesh_z_array)

	borders_mesh = [max_x, min_x, max_y, min_y, max_z, min_z]

	print('Maximum X:  ' + str(max_x) + '   km')
	print('Minimum X:  ' + str(min_x) + '   km')
	print('Maximum Y:  ' + str(max_y) + '   km')
	print('Minimum Y:  ' + str(min_y) + '   km')
	print('Maximum Z:  ' + str(max_z) + '   km')
	print('Minimum Z:  ' + str(min_z) + '   km')

	increment_in_x = np.abs(mesh_x_array[1] - mesh_x_array[0])

	increment_in_y = np.abs(np.unique(mesh_y_array)[1] - np.unique(mesh_y_array)[0])

	increment_in_z = np.abs(np.unique(mesh_z_array)[1] - np.unique(mesh_z_array)[0])

	x_mesh = np.arange(min_x, max_x + increment_in_x, increment_in_x)
	y_mesh = np.arange(min_y, max_y + increment_in_y, increment_in_y)
	z_mesh = np.arange(max_z, min_z, -increment_in_z)

	x_mesh_centers = x_mesh[:-1] + (increment_in_x / 2.0)
	y_mesh_centers = y_mesh[:-1] + (increment_in_y / 2.0)
	z_mesh_centers = z_mesh[:-1] - (increment_in_z / 2.0)

	mesh = np.meshgrid(x_mesh, y_mesh, z_mesh)
	mesh_center = np.meshgrid(x_mesh_centers, y_mesh_centers, z_mesh_centers)

	return mesh, mesh_center, x_mesh, y_mesh, z_mesh, x_mesh_centers, y_mesh_centers, z_mesh_centers, borders_mesh

def setup_material(material_data, material_names):

	"""A function to setup materials with the given material_data and material_names.
	It will automatically get rid of the air layers.
	
	Input:
	h5object: material_data - material h5 file dat ainput
	list: material_names - names of the material scraped from the py file.
					use the function read_uw_material_names_from_py_input to get this list.
	
	Output:
	array: material_array
	array: air_material_idx
	"""

	material_data_array = list(material_data['data'])

	material_array = []

	for i in range(0,len(material_data_array)):
		material_array.append(round(material_data_array[i][0]))

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

def setup_uw_data_array_PROJ_3D(data):

	"""A function to setup 3D underworld data array.
	"""

	data_list = list(data['data'])

	array = np.zeros(len(data_list))

	for i in range(0,len(data_list)):
		
		try:
			array[i] = data_list[i][0]
		except IndexError:
			array[i] = data_list[i]

	return array

def setup_uw_data_array_PROJ_2D(data):

	"""A function to setup 2D underworld data array.
	"""

	data_list = list(data['data'])

	array = np.zeros(len(data_list))

	for i in range(0,len(data_list)):

		array[i] = data_list[i][0]


	# array = np.array(list(data['data']))

	return array
