#!/usr/bin/env python3

import sys, csv
import numpy as np
import h5py

def read_h5_files(temp_h5, mesh_h5, material_h5, strain_h5):

	try:
		strain_file = h5py.File(strain_h5, 'r')
	except FileNotFoundError:
		print('The strain rate file is not found.')

	try:
		temp_file = h5py.File(temp_h5, 'r')

	except FileNotFoundError:
		print('The temperature file is not found.')
	try:
		mesh_file = h5py.File(mesh_h5, 'r')

	except FileNotFoundError:
		print('The mesh file is not found.')

	try:
		material_file = h5py.File(material_h5, 'r')

	except FileNotFoundError:
		print('The material file is not found.')
		
	
	return temp_file, mesh_file, material_file, strain_file
	
def read_uw_material_names_from_py_input(py_file):

	#Simple function for reading csv files and give out filtered output for given delimiter (delim)

	file_obj = open(py_file,'rt',encoding = "utf8") #Creating file object
	file_csv = csv.reader(file_obj) #Reading the file object with csv module, delimiter assigned to ','
	py_file_text = [] #Creating empty array to append data

	#Appending data from csb object
	for row in file_csv:
		py_file_text.append(row)

	#Filtering data for None elements read.
	for j in range(0,len(data)):
		py_file_text[j] = list(filter(None,data[j]))
	py_file_text = list(filter(None,data))

	
	idx_material_py = []

	for i in range(0,len(py_file_text)):
		if ('.add_material' in py_file_text[i][0]) == True:
			idx_material_py.append(i)
	
	material_names = []

	for i in idx_material_py:
		idx_local = py_file_text[i][0].find('=')
		material_names.append(py_file_text[i][0][:idx_local].strip())
	
	return material_names