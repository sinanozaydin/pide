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
	
	
def convert_3DModel_2_ModEM(file_out, conductivity_array, mesh, core_bounds = None, station_array = None, **kwargs):

	#This is a function that converts the constructed 3D geodynamic model into data points
	#that can be readable by the 3D MT modelling algorithm ModEM.
	"""
	Input Paramaters
	-----------------
	file_out: filename or full directory to output the result. Includes the format (e.g., csv)
	conductivity_array: calculated conductivity field in np.ndarray()
	mesh: mesh of the associated conductivity field.
	core_bounds: boundaries dictionary in  {'left','right','top','bottom'}	
	
	Keyword Arguments:
	core_mesh_size:
	num_horiz_bounds:
	horiz_bound_incr:
	num_vert_bounds:
	vert_bound_incr:
	"""
	
	from scipy.interpolate import griddata
	
	# core_mesh_size = kwargs.pop('core_mesh_size', mesh[0][0][1][0] - mesh[0][0][0][0])
	num_horiz_bounds = kwargs.pop('num_horiz_bounds', 8)
	horiz_bound_incr = kwargs.pop('horiz_bound_incr', 2)
	num_vert_bounds = kwargs.pop('num_vert_bounds', 9)
	vert_bound_incr = kwargs.pop('vert_bound_incr', 2)
	
	x_mesh = mesh[0]
	y_mesh = mesh[1]
	z_mesh = mesh[2]

	slice_len = len(x_mesh) * len(y_mesh)
	rho = np.zeros((len(z_mesh),slice_len))
	for i in range(0,len(z_mesh)):
		for j in range(0, slice_len):
			try:		
				rho[i][j] = conductivity_array[(i*slice_len)+j]
			except IndexError:
				raise IndexError('The mesh structure entered does not match the conductivity array. Be sure the entered format mesh = (x_mesh_centers,y_mesh_centers,z_mesh_centers) in tuples are correct.')
					
	#determining the indexes
	start_index_list = []
	sea_index_list = []
	for i in range(0,len(rho)):
		if np.all(rho[i] == rho[i][0]) == True:
			start_index_list.append(i)
		if (-999.0 in rho[i]) == False:
			sea_index_list.append(i)
	
	sea_index = sea_index_list[-1]
	
	for i in range(0,len(rho)):
		if i > sea_index:
			rho[i][(rho[i] == -999.0)] = 1e-14
	
	outs = np.float64(x_mesh[1] - x_mesh[0])
	xy_out = []
	for i in range(0, num_horiz_bounds):
		outs = outs * horiz_bound_incr
		xy_out.append(outs)
	xy_out = np.array(xy_out)
		
	outs_y = np.float64(y_mesh[1] - y_mesh[0])
	yx_out = []
	for i in range(0, num_horiz_bounds):
		outs_y = outs_y * horiz_bound_incr
		yx_out.append(outs_y)
	yx_out = np.array(yx_out)
	
	outs_z = np.float64(z_mesh[-2] - z_mesh[-1])
	z_out = []
	for i in range(0, num_vert_bounds):
		outs_z = outs_z * vert_bound_incr
		z_out.append(outs_z)
	outs_z = np.array(outs_z)
	
	x_core = np.ones(len(x_mesh)) * np.float64(x_mesh[1] - x_mesh[0])
	y_core = np.ones(len(y_mesh)) * np.float64(y_mesh[1] - y_mesh[0])
	z_core = np.ones(len(z_mesh)) * np.abs(np.float64(z_mesh[0] - z_mesh[1]))
	
	x_out = np.concatenate([xy_out[::-1],x_core,xy_out])
	y_out = np.concatenate([yx_out[::-1],y_core,yx_out])
	z_out = np.concatenate([z_core,z_out])
			
	rho_new = []
	
	for i in range(0,len(z_mesh)):
		rho_local = []
		for j in range(0,(num_horiz_bounds*((2*num_horiz_bounds)+len(x_core)))):
			rho_local.append(np.nan)
		for j in range(0,len(y_core)):
			for k in range(0,num_horiz_bounds):
				rho_local.append(np.nan)
			for k in range(j*len(x_core),j*len(x_core) + len(x_core)):
				rho_local.append(rho[i][k])
			for k in range(0,num_horiz_bounds):
				rho_local.append(np.nan)
		for j in range(0,(num_horiz_bounds*((2*num_horiz_bounds)+len(x_core)))):
			rho_local.append(np.nan)
		rho_new.append(np.array(rho_local))
	
	rho_new = np.array(rho_new)
	
	for i in range(0,num_vert_bounds):
		rho_new = np.insert(rho_new,0,np.ones(len(rho_new[-1])) * np.nan, axis = 0)
	
	for i in range(0,len(rho)):
		rho_new[i][0] = rho[i][0]
		rho_new[i][len(x_core) + 2*len(xy_out)] = rho[i][len(x_core)]
		rho_new[i][-(len(x_core) + 2*len(xy_out))] = rho[i][-len(x_core)]
		rho_new[i][-1] = rho[i][-1]
		
	for i in range(0,num_vert_bounds):
		rho_new[i] = rho_new[num_vert_bounds]
		
	rho_interp_array = rho_new.ravel()
	
	xi = np.cumsum(x_out)
	yi = np.cumsum(y_out)
	zi = np.cumsum(z_out)

	xi = np.insert(xi,0,0.0)
	yi = np.insert(yi,0,0.0)
	zi = np.insert(zi,0,0.0)
	
	xi_n = [((xi[i] - xi[i-1]) / 2.0) + xi[i-1] for i in range(1, len(xi))]
	yi_n = [((yi[i] - yi[i-1]) / 2.0) + yi[i-1] for i in range(1, len(yi))]
	zi_n = [((zi[i] - zi[i-1]) / 2.0) + zi[i-1] for i in range(1, len(zi))]
	
	zi_n = zi_n + z_mesh[-1] #adjusting the air layers
	
	zi_n = zi_n[::-1]
	
	x = []
	y = []
	z = []

	for i in range(0,len(zi_n)):
		for j in range(0,len(yi_n)):
			for k in range(0,len(xi_n)):
				x.append(xi_n[k])
				y.append(yi_n[j])
				z.append(zi_n[i])
				
	x = np.array(x)
	y = np.array(y)
	z = np.array(z)
	
	#getting the mask values for nan values to leave them out of the interpolation.
	mask = np.isnan(rho_interp_array)
	mesh_out = np.meshgrid(xi_n,yi_n,zi_n)
	
	#interpolating the data using griddata
	rho__ = griddata((x[~mask],y[~mask],z[~mask]), rho_interp_array[~mask], (mesh_out[0],mesh_out[1],mesh_out[2]), method = 'nearest')
	
	#rearranging the griddata output to put in a .rho file more easily.
	rho_write = []
	for i in range(0,len(zi_n)):
		local_rho = []
		for j in range(0,len(yi_n)):
			for k in range(0,len(xi_n)):
				local_rho.append(rho__[k][j][i])
		rho_write.append(local_rho)
				
	rho_write = np.array(rho_write)
	
	#reversing course for writing files
	rho_write = rho_write[::-1]
	
	#finding the index where the layer is all air
	s_idx = []
	for i in range(0,len(rho_write)):
		if np.all(rho_write[i] == rho_write[i][0]) == True:
			s_idx.append(i)
			
	#Getting rid of the layers only contains air
	rho_write = rho_write[s_idx[-1]+1:]
	z_out = z_out[s_idx[-1]+1:]
		
	rho_write = np.log(1.0 / rho_write) #ModEM rho format with natural logarithm of resistivityrho_
	
	lines = ["#  3D MT model written by ModEM in WS format by pide\n"]
	lines.append('   '+str(len(x_out))+ '   ' +str(len(y_out))+ '   '+str(len(z_out))+ ' 0    LOGE\n')
	
	line = np.array2string(x_out*1e3, separator=' ', max_line_width=np.inf, formatter={'all': lambda x: f'{x:10.3f}'})
	lines.append('  ' + line[1:-1] + '\n')
	
	line = np.array2string(y_out*1e3, separator=' ', max_line_width=np.inf, formatter={'all': lambda x: f'{x:10.3f}'})
	lines.append('  ' + line[1:-1] + '\n')
	
	line = np.array2string(z_out*1e3, separator=' ', max_line_width=np.inf, formatter={'all': lambda x: f'{x:10.3f}'})
	lines.append('  ' + line[1:-1] + '\n')
	lines.append('\n')

	for i in range(0,len(rho_write)):
		for j in range(0,len(rho_write[i]),len(x_out)):
			line = np.array2string(rho_write[i][j:(j + len(x_out))], separator='  ', max_line_width=np.inf, formatter={'all': lambda x: f'{x:.5E}'})
			lines.append('  ' + line[1:-1] + '\n')
		lines.append('\n')
		
	center_z = 0.0
	center_east = -0.5 * (np.sum(np.abs(x_out)))
	center_north = -0.5 * (np.sum(np.abs(y_out)))
	grid_center = ['    ',str(center_north),'   ',str(center_east),'   ',str(center_z),'\n']
	grid_center = ''.join(grid_center)
	lines.append(grid_center)
	lines.append('    0.000\n')
	
	filesave_composition = open(file_out ,'w')
	filesave_composition.writelines(lines)
	filesave_composition.close()
	print('File for the converted model is written as: ' + file_out)
	#Finding the uppermost layer with no nan values
	
def create_ModEM_fwd_file(file_out, input_rho, station_location_arrays = None, freq_start=100, freq_end=1e-4, num_freq=120):

	try:
		file_name = input_rho[input_rho.rfind('/')+1:-3]
	except:
		file_name = input_rho[:-3]
		
	from pide.mt.mt_model_read import read_ModEM_rho
	from pide.utils.gis_tools import utm_to_lat_lon, lat_lon_to_utm, get_utm_zone_number
	import decimal

	rho, mesh_centers_x_array, mesh_centers_y_array, z_mesh_center = read_ModEM_rho(rho_file_path=input_rho)
	
	#finding the unique grid locations
	mesh_x_uniq = np.unique(mesh_centers_x_array)
	mesh_y_uniq = np.unique(mesh_centers_y_array)
	
	x_grid = [mesh_x_uniq[i] - mesh_x_uniq[i-1] for i in range(1,len(mesh_x_uniq))]
	y_grid = [mesh_y_uniq[i] - mesh_y_uniq[i-1] for i in range(1,len(mesh_y_uniq))]
	
	#finding the core start and end
	idx_x_count = 0
	for i in range(1,len(x_grid)):
		if x_grid[i] == x_grid[i-1]:
			idx_x_count = idx_x_count + 1
			if idx_x_count == 1:
				idx_x_start = i-1
			else:
				idx_x_end = i

	idx_y_count = 0
	for i in range(1,len(y_grid)):
		if y_grid[i] == y_grid[i-1]:
			idx_y_count = idx_y_count + 1
			if idx_y_count == 1:
				idx_y_start = i-1
			else:
				idx_y_end = i
	
	if station_location_arrays is None:
		raise KeyError('Station locations has to be entered to create a file')
	else:
		if np.mean([len(station_location_arrays[0]),len(station_location_arrays[1]),len(station_location_arrays[2])]) == len(station_location_arrays[0]):
			station_location_arrays = np.array(station_location_arrays) * 1e3
		else:
			raise IndexError('The length of the station_locations_array are not the same.')
			
	utm_zone = get_utm_zone_number(longitude = 0.0)
	mc_x, mc_y = lat_lon_to_utm(0,0,zone_number=utm_zone)
	x_loc = station_location_arrays[0] + mc_x
	y_loc = station_location_arrays[1] + mc_y
	lat_degrees_list, lon_degrees_list = utm_to_lat_lon(x_loc,y_loc,zone_number=utm_zone)
	
	#creating frequency array
	freq_array = np.logspace(np.log10(freq_start), np.log10(freq_end),num_freq)
		
	dat_lines = []
	dat_lines.append('# ModEM impedance responses' +  '\n')
	dat_lines.append('# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n')
	dat_lines.append('> Full_Impedance\n')
	dat_lines.append('> exp(+i\omega t)\n')
	dat_lines.append('> [mV/km]/[nT]\n')
	dat_lines.append('> 0.00\n')
	dat_lines.append('> ' + str(0.00) + ' ' + str(0.00) + '\n' )
	dat_lines.append('> ' + str(len(freq_array)) + ' ' + str(len(station_location_arrays[0])) + '\n')
	for i in range(0,len(station_location_arrays[0])):
		for j in range(0,len(freq_array)):
			for k in range(0,4):
				if k == 0:
					dum = 'ZXX'
					dum2 = '% .5E' % decimal.Decimal(str(0.0))
					dum3 = '% .5E' % decimal.Decimal(str(0.0))
					dum4 = '% .5E' % decimal.Decimal(str(0.0))
				elif k == 1:
					dum = 'ZXY'
					dum2 = '% .5E' % decimal.Decimal(str(0.0))
					dum3 = '% .5E' % decimal.Decimal(str(0.0))
					dum4 = '% .5E' % decimal.Decimal(str(0.0))

				elif k == 2:
					dum = 'ZYX'
					dum2 = '% .5E' % decimal.Decimal(str(0.0))
					dum3 = '% .5E' % decimal.Decimal(str(0.0))
					dum4 = '% .5E' % decimal.Decimal(str(0.0))

				elif k == 3:
					dum = 'ZYY'
					dum2 = '% .5E' % decimal.Decimal(str(0.0))
					dum3 = '% .5E' % decimal.Decimal(str(0.0))
					dum4 = '% .5E' % decimal.Decimal(str(0.0))


				dat_lines.append(str('%.5E' % decimal.Decimal(1.0/freq_array[j])) + '  ' + ('%-15s' % ('C_' + str(i))) + '  ' + str('% 7.4f' % lat_degrees_list[i]) +\
					'  ' + str('% 7.4f' % lon_degrees_list[i]) + '  ' + str('% 12.3f' % (station_location_arrays[0][i])) + '  ' +\
					str('% 12.3f' %(station_location_arrays[1][i])) + '  ' + str('%8.3f' % float(station_location_arrays[2][i])) + '  ' +\
					 dum + '  ' + dum2 + '  ' + dum3 + '  ' + dum4 + '\n')

	dat_lines.append('#\n')

	filesave = open(file_out,'w')
	filesave.writelines(dat_lines)
	filesave.close
	print('The data (.dat) file is written at : ' + os.getcwd())
	
def get_station_elevation_ModEM_rho(input_rho, station_location_arrays = None, air_conductivity = 1e-14):

	from pide.mt.mt_model_read import read_ModEM_rho
		
	if station_location_arrays is None:
		raise KeyError('Station locations has to be entered to create a file')
	else:
		if np.mean([len(station_location_arrays[0]),len(station_location_arrays[1])]) == len(station_location_arrays[0]):
			station_location_arrays = np.array(station_location_arrays) * 1e3
		else:
			raise IndexError('The length of the station_locations_array are not the same.')
	
	#reading the rho file.
	rho, x_mesh_center, y_mesh_center, z_mesh_center = read_ModEM_rho(rho_file_path=input_rho)

	#finding the unique grid locations
	mesh_x_uniq = np.unique(x_mesh_center)
	mesh_y_uniq = np.unique(y_mesh_center)
	mesh_z_uniq = np.unique(z_mesh_center)

	z_grid = [mesh_z_uniq[i] - mesh_z_uniq[i-1] for i in range(1,len(mesh_z_uniq))]
	z_depth = np.array(mesh_z_uniq[:-1]) - (np.array(z_grid)/2.0)
	
	air_res = np.log10(1.0/air_conductivity)

	depth_sol = []

	for i in range(0,len(station_location_arrays[0])):
	
		idx_x = (np.abs(np.asarray(mesh_x_uniq) - station_location_arrays[0][i])).argmin()
		idx_y = (np.abs(np.asarray(mesh_y_uniq) - station_location_arrays[1][i])).argmin()

		rho_profile = np.log10([rho[j][idx_y][idx_x] for j in range(0,len(rho))])

		idx_z_list = [j for j in range(0,len(rho_profile)) if (round(rho_profile[j],4) == air_res)]
		
		try:
			idx_z = idx_z_list[-1] + 1
		except IndexError:
			idx_z = 0

		depth_sol.append(z_depth[idx_z])
			
	return (station_location_arrays[0]*1e-3,station_location_arrays[1]*1e-3,np.array(depth_sol)*1e-3)
		

	
	
	
	
	
	
	
	