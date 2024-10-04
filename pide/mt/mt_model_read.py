import numpy as np
from pide.utils.utils import read_csv

def read_ModEM_rho(rho_file_path):

	"""A functions that reads ModEM file rho format and gives relevant data
	
	Input:
	str: rho_file_path
	
	Output:
	array: rho - 3D matrix of resistivity values || in ohm.m
	array: mesh_centers_x_array - locations of the mesh centers in x direction for each cell
	array: mesh_centers_y_array - locations of the mesh centers in y direction for each cell
	array: z_mesh_center - locations of the mesh centers in z direction for each cell
	"""

	ModEM_rho_data = read_csv(filename = rho_file_path, delim = ' ')

	x_num = int(ModEM_rho_data[1][0])
	y_num = int(ModEM_rho_data[1][1])
	z_num = int(ModEM_rho_data[1][2])

	x_grid = np.array(ModEM_rho_data[2],dtype=float)
	y_grid = np.array(ModEM_rho_data[3],dtype=float)
	z_grid = np.array(ModEM_rho_data[4],dtype=float)


	lenxgrid = len(x_grid)
	lenygrid = len(y_grid)
	lenzgrid = len(z_grid)
	
	
	rho = []

	for k in range(5,len(ModEM_rho_data) - 2 ,y_num):
		rhoy = []
		for z in range(k, k + y_num):
			rhox = []
			for l in range(0,x_num):
				rhox.append(float(ModEM_rho_data[z][l]))
			rhoy.append(rhox)

		rho.append(rhoy)

	rho = np.exp(np.array(rho))

	
	z_depth = np.array([0.0])
	z_grid = np.cumsum(z_grid)
	
	z_depth = np.append(z_depth,z_grid)
	z_mesh_center = []
	for i in range(1,len(z_depth)):
		z_mesh_center.append((((z_depth[i] - z_depth[i-1]) / 2.0) + z_depth[i-1]))
	z_mesh_center = np.array(z_mesh_center)

	x_grid_cum = []
	y_grid_cum = []
	mid_point_x = int(len(x_grid) / 2.0)
	mid_point_y = int(len(y_grid) / 2.0)
	
	if len(x_grid) %2 == 0:
		mid_point_x = int(len(x_grid) / 2.0)
		beg_x = np.sum(x_grid[:mid_point_x]) * -1

	elif len(x_grid) %2 != 0:
		mid_point_x = int(len(x_grid) / 2.0) + 1
		beg_x = np.sum(x_grid[:mid_point_x]) * -1 + (x_grid[mid_point_x] / 2.0)

	if len(y_grid) %2 == 0:
		mid_point_y = int(len(y_grid) / 2.0)
		beg_y = np.sum(y_grid[:mid_point_y]) * -1
	elif len(y_grid) %2 != 0:
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

	x_grid_cum = x_grid_cum[::-1]
	x_grid = x_grid[::-1]

	#Creating x and y mesh centers to find the profile locations

	x_mesh_center = []
	y_mesh_center = []

	for i in range(1,len(x_grid_cum)):
		x_mesh_center.append(((x_grid_cum[i] - x_grid_cum[i-1]) / 2.0) + x_grid_cum[i-1])
	for i in range(1,len(y_grid_cum)):
		y_mesh_center.append(((y_grid_cum[i] - y_grid_cum[i-1]) / 2.0) + y_grid_cum[i-1])

	mesh_centers = np.meshgrid(x_mesh_center,y_mesh_center)
	
	mesh_centers_x_array = mesh_centers[1].flatten()
	mesh_centers_y_array = mesh_centers[0].flatten()
	
	
	return rho, mesh_centers_x_array, mesh_centers_y_array, z_mesh_center
