
import os

from pide.mt.mt_model_read import read_ModEM_rho
from pide.utils.utils import read_csv, associate_coordinates

import numpy as np
from mtpy.utils import gis_tools as gis_tools
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

import pide

rho_file_path = "/home/sinan/src/SEEL/examples/scripts/e_aust/J1_mean_model.rho"

rho, mesh_x, mesh_y, mesh_z = read_ModEM_rho(rho_file_path=rho_file_path)

#mc_lat and lon from the dat file
mc_lat = -34.705091
mc_lon = 145.677690

#converting mc_lat and lon to x y as to get reference for the utm coordinates
utm_no = gis_tools.get_epsg(mc_lat,mc_lon)
stuff_center = gis_tools.project_point_ll2utm(mc_lat, mc_lon, utm_zone = 34, epsg = utm_no)

mc_x = stuff_center[0]
mc_y = stuff_center[1]

#creating utm coordinate x and y from mesh_x and mesh_y
x = mesh_x + mc_x
y = mesh_y + mc_y


#converting utm corrdinates back to latitude longitude
stuff = gis_tools.project_point_utm2ll(x, y, utm_zone = 34, epsg = utm_no)
		
lon_modem =  np.zeros(len(stuff))
lat_modem = np.zeros(len(stuff))

for i in range(0,len(stuff)):
	lat_modem[i] = stuff[i][0]
	lon_modem[i] = stuff[i][1]

	
#finding 50 km slice from rho
idx_50 = (np.abs((mesh_z/1e3) - 50)).argmin()
depth_str = str(mesh_z[idx_50]/1e3) + ' km'

rho_50 = rho[idx_50].flatten()

#getting the temperature file from Manassero et al. (2024)
filename_T = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'e_aust', 'temp_Manassero2024.dat')
data_T = read_csv(filename=filename_T, delim = ' ')
data_T = np.array(data_T, dtype = float)

data_T_depth = data_T[:,3]
data_T_temp = data_T[:,4] + 273.15 #converting it to K
data_T_lat = data_T[:,1]
data_T_lon = data_T[:,2]




#interpolating temperature at 48.7 kilometer dept
T_40 = data_T_temp[np.where(data_T_depth==-40)[0]]
T_60 = data_T_temp[np.where(data_T_depth==-60)[0]]

data_lat = data_T_lat[np.where(data_T_depth==-40)[0]]
data_lon = data_T_lon[np.where(data_T_depth==-40)[0]]

#converting mc_lat and lon to x y as to get reference for the utm coordinates
stuff = gis_tools.project_point_ll2utm(data_lat, data_lon, utm_zone = 34, epsg = utm_no)

x_T =  np.zeros(len(stuff))
y_T = np.zeros(len(stuff))

for i in range(0,len(stuff)):
	x_T[i] = stuff[i][0]
	y_T[i] = stuff[i][1]

T_48 = np.zeros(len(T_40))

for i in range(0,len(T_40)):

	T_48[i] = np.interp(48.7, [40,60], [T_40[i], T_60[i]])


#2D interpolation of Temperature field to make 

x_new = np.arange(np.amin(x_T),np.amax(x_T),5e3) #every 5_000 m
y_new = np.arange(np.amin(y_T),np.amax(y_T),5e3) #every 5_000 m

coords_mesh = np.meshgrid(x_new, y_new)

points_interp = np.column_stack((x_T, y_T))

T_48_interp = griddata(points_interp, T_48, (coords_mesh[0], coords_mesh[1]), method = 'cubic')

T_48_interp = T_48_interp.flatten()
coords_x = coords_mesh[0].flatten()
coords_y = coords_mesh[1].flatten()

#Applying mask for nan values at T due to interpolation out of bounds
mask = np.isnan(T_48_interp)
T_48_i = T_48_interp[~mask]
coords_x = coords_x[~mask]
coords_y = coords_y[~mask]

#associateding coordinates for with

idx_ = associate_coordinates(sample_x = coords_x, sample_y = coords_y, target_x=x, target_y=y, num_cpu = 5)

p_obj = pide.pide()
p_obj.set_temperature(T_48_i[idx_])
p_obj.set_pressure(48.7 / 33.0)
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.25,cpx = 0.1, garnet= 0.05)

p_obj.list_mineral_econd_models('cpx')

p_obj.set_mineral_conductivity_choice(ol = 4, opx = 0, garnet = 0, cpx = 6)
p_obj.set_solid_phs_mix_method(1) #HS lower bound

cond = p_obj.calculate_conductivity()


logcond = np.log10(cond)
logrho = np.log10(1.0/rho_50)

diff_log = logrho-logcond

#Writing the results
lines = ['Latitude,Longitude,DiffLogCond,SyntheticLogCond,ModelLogCond\n']

for i in range(len(diff_log)):
	
	lines.append(','.join((str(lat_modem[i]),str(lon_modem[i]),str(diff_log[i]),str(logcond[i]),str(logrho[i]) + '\n')))
	
with open('LogCond.csv','w') as file:
	file.writelines(lines)

"""
fig = plt.figure()
ax = plt.subplot(111)

cax = ax.scatter(x, y, c = diff_log, cmap = 'Spectral_r', marker = 'o')
cax.set_clim(np.amin(diff_log), np.amax(diff_log))
cbar_xx = fig.colorbar(cax,boundaries=np.linspace(np.amin(diff_log), np.amax(diff_log)),orientation="vertical", pad=0.05,
		 ax = ax)
ax.plot()
plt.show()


"""