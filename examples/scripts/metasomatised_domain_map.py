
import os

from pide.mt.mt_model_read import read_ModEM_rho
from pide.utils.utils import read_csv

import numpy as np
from mtpy.utils import gis_tools as gis_tools
import matplotlib.pyplot as plt

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
data_T_lon = data_T[:2]

T_40 = data_T_temp[np.where(data_T_depth==-40)[0]]
T_60 = data_T_temp[np.where(data_T_depth==-60)[0]]

T_48 = np.zeros(len(T_40))

for i in range(0,len(T_40)):

	T_48[i] = np.interp(48.7, [40,60], [T_40[i], T_60[i]])




fig = plt.figure()
ax = plt.subplot(111)

cax = ax.scatter(lon_modem, lat_modem, c = np.log10(rho[16].flatten()), cmap = 'Spectral', marker = 'o')
cax.set_clim(0,4)
cbar_xx = fig.colorbar(cax,boundaries=np.linspace(0,4),orientation="vertical", pad=0.05,
		 ticks = [0,1,2,3,4], ax = ax)
ax.plot()
plt.show()


