import numpy as np
from pide.utils.utils import read_csv

def read_ModEM_rho(rho_file_path):

	"""
	A function that reads the ModEM format rho file.
	
	__author__ = "Sinan Ozaydin | sinan.ozaydin@protonmail.com"

	Parameters
	----------
	rho_file_path: str
		Path to a single rho file.
		
	Returns
	-------
	rho : ndarray of float
		Apparent resistivity in ohm-meters.
	
	rho_dict : dict
		A dictionary containing mesh boundary values and other positional parameters with the following keys:
		
		- 'min_x_grid_cum' : float  
			Minimum value of the cumulative X grid.
	
		- 'max_x_grid_cum' : float  
			Maximum value of the cumulative X grid.
	
		- 'min_y_grid_cum' : float  
			Minimum value of the cumulative Y grid.
	
		- 'max_y_grid_cum' : float  
			Maximum value of the cumulative Y grid.
	
		- 'x_mesh_center' : float  
			Center coordinate of the mesh in the X direction.
	
		- 'y_mesh_center' : float  
			Center coordinate of the mesh in the Y direction.
	
		- 'z_mesh_center' : float  
			Center coordinate of the mesh in the Z direction.
	
		- 'mesh_centers_x_array' : np.ndarray  
			Array of X-direction mesh center coordinates.
	
		- 'mesh_centers_y_array' : np.ndarray  
			Array of Y-direction mesh center coordinates.
	
	"""

	ModEM_rho_data = read_csv(filename = rho_file_path, delim = ' ')

	x_num = int(ModEM_rho_data[1][0])
	y_num = int(ModEM_rho_data[1][1])
	z_num = int(ModEM_rho_data[1][2])
	
	x_grid = np.array([float(x) for x in ModEM_rho_data[2]])
	y_grid = np.array([float(x) for x in ModEM_rho_data[3]])
	z_grid = np.array([float(x) for x in ModEM_rho_data[4]])
	
	lenxgrid = len(x_grid)
	lenygrid = len(y_grid)
	lenzgrid = len(z_grid)

	rho = []
	rho_chain = []

	for k in range(5,len(ModEM_rho_data) - 2 ,y_num):
		rhoy = []

		for z in range(k, k + y_num):
			rhox = []
			for l in range(x_num-1,-1,-1):
				rhox.append(float(ModEM_rho_data[z][l]))
				rho_chain.append(float(ModEM_rho_data[z][l]))
			rhoy.append(rhox)

		rho.append(rhoy)

	rho = np.exp(np.asarray(rho))
	rho_chain = np.exp(np.array(rho_chain))
	
	z_depth = np.array([0.0])
	z_grid = np.cumsum(z_grid)
	z_depth = np.append(z_depth,z_grid)
	z_mesh_center = []
	for i in range(1,len(z_depth)):
		z_mesh_center.append((((z_depth[i] - z_depth[i-1]) / 2.0) + z_depth[i-1]) / 1000.0)

	z_mesh_center = np.array(z_mesh_center)
	
	x_grid_cum = []
	y_grid_cum = []
	mid_point_x = int(len(x_grid) / 2.0)
	mid_point_y = int(len(y_grid) / 2.0)

	if len(x_grid) % 2 == 0:
		mid_point_x = int(len(x_grid) / 2.0)
		beg_x = np.sum(x_grid[:mid_point_x]) * -1
	elif len(x_grid) % 2 != 0:
		mid_point_x = int(len(x_grid) / 2.0) + 1
		beg_x = np.sum(x_grid[:mid_point_x]) * -1 + (x_grid[mid_point_x] / 2.0)

	if len(y_grid) % 2 == 0:
		mid_point_y = int(len(y_grid) / 2.0)
		beg_y = np.sum(y_grid[:mid_point_y]) * -1
	elif len(y_grid) % 2 != 0:
		mid_point_y = int(len(y_grid) / 2.0) + 1
		beg_y = np.sum(y_grid[:mid_point_y]) * -1 + (y_grid[mid_point_y] / 2.0)
	
	x_grid_cum.append(beg_x)
	for i in range(0,lenxgrid):
		beg_x += x_grid[i]
		x_grid_cum.append(beg_x)
	y_grid_cum.append(beg_y)
	for i in range(0,lenygrid):
		beg_y += y_grid[i]
		y_grid_cum.append(beg_y)

	min_x_grid_cum = np.amin(x_grid_cum)
	max_x_grid_cum = np.amax(x_grid_cum)
	min_y_grid_cum = np.amin(y_grid_cum)
	max_y_grid_cum = np.amax(y_grid_cum)
	
	#Creating x and y mesh centers to find the profile locations

	x_mesh_center = []
	y_mesh_center = []

	for i in range(1,len(x_grid_cum)):
		x_mesh_center.append(((x_grid_cum[i] - x_grid_cum[i-1]) / 2.0) + x_grid_cum[i-1])
	for i in range(1,len(y_grid_cum)):
		y_mesh_center.append(((y_grid_cum[i] - y_grid_cum[i-1]) / 2.0) + y_grid_cum[i-1])

	x_mesh_center = np.array(x_mesh_center) * 1e-3
	y_mesh_center = np.array(y_mesh_center) * 1e-3
	
	mesh_centers = np.meshgrid(x_mesh_center,y_mesh_center)
	
	mesh_centers_x_array = mesh_centers[1].flatten()
	mesh_centers_y_array = mesh_centers[0].flatten()
	
	return rho, {'min_x_grid_cum': min_x_grid_cum, 'max_x_grid_cum': max_x_grid_cum, 'min_y_grid_cum': min_y_grid_cum, 'max_y_grid_cum': max_y_grid_cum,
	'x_mesh_center': x_mesh_center, 'y_mesh_center': y_mesh_center, 'z_mesh_center': z_mesh_center,
	'mesh_centers_x_array': mesh_centers_x_array, 'mesh_centers_y_array': mesh_centers_y_array,
	'x_grid_cum': x_grid_cum, 'y_grid_cum': y_grid_cum, 'x_grid':x_grid, 'y_grid': y_grid} 

def read_ModEM_dat(dat_file, reading_for = 'data'):

	#Reading ModEM dat file to get coordinates of stations and model centre.

	ModEM_dat_data = read_csv(dat_file,delim = ' ')

	station_lat = []
	station_lon = []
	station_name = []
	station_posx = []
	station_posy = []

	dash_found = False
	
	for rows in range(0,15):
		if ModEM_dat_data[rows][0] == '>':
			start_idx = rows
	start_idx = start_idx + 1
	
	mc_lat = float(ModEM_dat_data[start_idx-2][1])
	mc_lon = float(ModEM_dat_data[start_idx-2][2])
	
	for row in range(start_idx,len(ModEM_dat_data)):
		if ModEM_dat_data[row][0] == '#' :
			limitlines = row-1
			dash_found = True

	if dash_found == False:
		limitlines = len(ModEM_dat_data)

	for row in range(start_idx,limitlines):
		if row == start_idx:
			station_name.append(ModEM_dat_data[row][1])
			station_lat.append(float(ModEM_dat_data[row][2]))
			station_lon.append(float(ModEM_dat_data[row][3]))
			station_posx.append(float(ModEM_dat_data[row][5]))
			station_posy.append(float(ModEM_dat_data[row][4]))
		else:
			if (ModEM_dat_data[row][4] != ModEM_dat_data[row-1][4]) or (ModEM_dat_data[row][5] != ModEM_dat_data[row-1][5]):
				station_name.append(ModEM_dat_data[row][1])
				station_lat.append(float(ModEM_dat_data[row][2]))
				station_lon.append(float(ModEM_dat_data[row][3]))
				station_posx.append(float(ModEM_dat_data[row][5]))
				station_posy.append(float(ModEM_dat_data[row][4]))

	st_item_list = []
	
	for i in range(0,len(station_name)):
	
			st_item_list.append(str(i) + ' - ' + station_name[i])
	
	station_lat = np.array(station_lat)
	station_lon = np.array(station_lon)
	station_posx = np.array(station_posx)
	station_posy = np.array(station_posy)
	
	if reading_for == 'positioning':
		
		return {'station_lat': station_lat, 'station_lon':station_lon, 'station_posx': station_posx,
		'station_posy':station_posy, 'mc_lat': mc_lat, 'mc_lon': mc_lon, 'station_name': station_name}
