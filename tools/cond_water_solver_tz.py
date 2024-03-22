import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from functools import partial
from scipy.io import savemat, loadmat

import os,sys,csv
import netCDF4

def read_csv(filename,delim, filter_rows =True):

	#Simple function for reading csv files and give out filtered output for given delimiter (delim)

	file_obj = open(filename,'rt',encoding = "utf8") #Creating file object
	file_csv = csv.reader(file_obj,delimiter = delim) #Reading the file object with csv module, delimiter assigned to ','
	data = [] #Creating empty array to append data

	#Appending data from csb object
	for row in file_csv:
		data.append(row)

	#Filtering data for None elements read.
	for j in range(0,len(data)):
		if filter_rows == True:
			data[j] = list(filter(None,data[j]))
	data = list(filter(None,data))

	return data
	
def write_slice(p,t,water,cond, lat, lon, depth, filename):

	lines = ['Latitude, Longitude,Depth[km],P[GPa],T[K],Water[ppm], cond[Sm]\n']
	
	for i in range(0,len(p)):
		
		if lon[i] > 180.0:
			lines.append(','.join((str(lat[i]), str(lon[i] - 360.0), str(depth), str(p[i]),str(t[i]),str(water[i]),str(cond[i]) + '\n')))
		else:
			lines.append(','.join((str(lat[i]), str(lon[i]), str(depth), str(p[i]),str(t[i]),str(water[i]),str(cond[i]) + '\n')))
		
	filesave = open(filename,'w')
	filesave.writelines(lines)
	filesave.close()
	
	print('Slice file is written at : ' + filename)
	
core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide')

sys.path.append(core_path_ext)
import pide
from inversion import conductivity_solver_single_param


filename_global_em = "/home/sinan/Desktop/Research/Transition_Zone_Ben/GlobalEM-2015-02x02.nc"

em_nc = netCDF4.Dataset(filename_global_em, 'r')

lon_model = np.asarray(em_nc.variables['longitude'][:])
lat_model = np.asarray(em_nc.variables['latitude'][:])
depth_model = np.asarray(em_nc.variables['depth'][:])
sigma = np.array(em_nc.variables['sigma'][:])

index_tz_start = 21
index_tz_end = 31

filename_phase_file = "/home/sinan/Desktop/Research/Transition_Zone_Ben/Phase_proportion.csv"
data_phases = read_csv(filename = filename_phase_file, delim = ',')

depth_of_interest = depth_model[index_tz_start:index_tz_end]

p = []
depth_p = []
t_p = []
cpx = []
garnet = []
wad = []
ring = []
pv = []

for i in range(1,len(data_phases)):

	if (float(data_phases[i][1]) > 410.0) and (float(data_phases[i][1]) <= 660.0):
		
		if any(element < 1.5 for element in abs(float(data_phases[i][1]) - depth_of_interest)):
			p.append(float(data_phases[i][0]) / 1e4) #converting to gpa
			depth_p.append(float(data_phases[i][1]))
			t_p.append(float(data_phases[i][2])) # in kelvin
			cpx.append(float(data_phases[i][13]) + float(data_phases[i][14]) + float(data_phases[i][15]) + float(data_phases[i][17]))
			garnet.append(float(data_phases[i][18]) + float(data_phases[i][19]) + float(data_phases[i][20]) + float(data_phases[i][21]) + float(data_phases[i][22]))
			wad.append(float(data_phases[i][26]) + float(data_phases[i][27]))
			ring.append(float(data_phases[i][28]) + float(data_phases[i][29]))
			pv.append(float(data_phases[i][23]))
			
p = np.array(p)
depth_p = np.array(depth_p)
t_p = np.array(t_p)
cpx = np.array(cpx)
garnet = np.array(garnet)
wad = np.array(wad)
ring = np.array(ring)
pv = np.array(pv)

#reading katsura mantle potential temperature gradient

filename_temp_grad_file = "/home/sinan/Desktop/Research/Transition_Zone_Ben/temperature_models/katsura_2020_grads.csv"
temp_grad_data = read_csv(filename = filename_temp_grad_file, delim = ',')

katsura_depth = []
katsura_grad = []

for i in range(1,len(temp_grad_data)):

	katsura_depth.append(float(temp_grad_data[i][0]))
	katsura_grad.append(float(temp_grad_data[i][1]))
	
#reading waszek global mantle potential temperature diferences.

temp_model_pot = "/home/sinan/Desktop/Research/Transition_Zone_Ben/temperature_models/410660/thermalmodel.txt"
temp_model_pot_data = read_csv(filename = temp_model_pot, delim = ' ', filter_rows=False)

waszek_t = []
waszek_t_list = []
waszek_lat = []
waszek_lon = []

waszek_depth = list(katsura_depth)

waszek_depth.insert(0, 0.0) #general depth
waszek_depth = np.array(waszek_depth)

#combining the waszek potential temperature map and Katsura2022 average adiabatic gradients.
for i in range(1,len(temp_model_pot_data)):

	waszek_lat.append(float(temp_model_pot_data[i][0]))
	waszek_lon.append(float(temp_model_pot_data[i][1]))
	waszek_t.append(float(temp_model_pot_data[i][2]))
	t_local = [float(temp_model_pot_data[i][2])]
	for j in range(0,len(katsura_grad)):
		
		if j == 0:
			t_local.append((katsura_depth[j] * katsura_grad[j]) + t_local[0])
			
		else:
			t_local.append(((katsura_depth[j] - katsura_depth[j-1]) * katsura_grad[j]) + t_local[j])
			
	waszek_t_list.append(np.array(t_local))
	
waszek_t_list = np.array(waszek_t_list)
waszek_lat = np.array(waszek_lat)
waszek_lon = np.array(waszek_lon)

#finding the corresponding depth indexes between depth_p and waszek_depth

idx_depth_list = []
for i in range(0,len(depth_p)):
	idx_depth_ = (np.abs(waszek_depth-depth_p[i])).argmin()
	if i != 0:
		if (idx_depth_ == idx_depth_list[-1]):
			idx_depth_ = idx_depth_ + 1
	idx_depth_list.append(idx_depth_)


#getting rid of the little deviation from 1 and subtracting from the major silicate phase
tot = cpx + garnet + wad + ring + pv
deficit = tot - 1.0

for i in range(0,len(tot)):

	if wad[i] == 0.0:
		ring[i] = ring[i] - deficit[i]
	else:
		wad[i] = wad[i] - deficit[i]
		
tot = cpx + garnet + wad + ring + pv

#setting up conductivity arrays
cond_arrays = []
lat_arrays = []
lon_arrays = []

for i in range(index_tz_start, index_tz_end):

	cond_slice = []
	for j in range(0,len(sigma[i])):
		for k in range(0,len(sigma[i][j])):
			cond_slice.append(10**sigma[i][j][k])
			
			if i == index_tz_start:
				lat_arrays.append(lat_model[j])
				if lon_model[k] > 180.0:
					lon_arrays.append(lon_model[k] - 360.0)
				else:
					lon_arrays.append(lon_model[k])
				
	cond_arrays.append(np.array(cond_slice))
	
cond_arrays = np.array(cond_arrays)
lat_arrays = np.array(lat_arrays)
lon_arrays = np.array(lon_arrays)

"""
#determining the corresponding coordinate indexes of temperature field generated with the resistivity matrix.

def _associate_coordinates_(index, lat_arrays, lon_arrays, waszek_lat, waszek_lon):

	idx_lat_target = (np.abs(waszek_lat-lat_arrays[index])).argmin()
	idx_lat_lists = [idx for idx, value in enumerate(waszek_lat) if value == waszek_lat[idx_lat_target]]
	idx_ = (np.abs(waszek_lon[idx_lat_lists]-lon_arrays[index])).argmin()
	idx_final = idx_lat_lists[idx_]
	
	return idx_final

#parallelisation to save time
num_cpu = 6

index_list = range(0,len(lat_arrays))

with multiprocessing.Pool(processes=num_cpu) as pool:
	
	process_item_partial = partial(_associate_coordinates_, lat_arrays = lat_arrays, lon_arrays = lon_arrays, waszek_lat = waszek_lat,
	waszek_lon = waszek_lon)
	
	c = pool.map(process_item_partial, index_list)
				
idx_temp_array = [x for x in c]
idx_temp_array = np.array(idx_temp_array)

T_global = waszek_t_list[idx_temp_array]
savemat('T_global_waszek.mat', {'matrix': T_global})

fig = plt.figure()
ax = plt.subplot(111)
ax.scatter(lon_arrays,lat_arrays,c = T_global[:,0], cmap = 'coolwarm', marker = 's',
		 linewidth = 0.05, edgecolor = 'k')
plt.show()

"""
#loading the temperature array created in the commented section above in the previous run
mat_contents = loadmat('T_global_waszek.mat')
T_global = mat_contents['matrix']

#creating pide object
pide_obj = pide.pide()
pide_obj.list_mineral_econd_models('garnet')
pide_obj.list_mineral_econd_models('cpx')
pide_obj.list_mineral_econd_models('perov')
pide_obj.list_mineral_econd_models('rwd_wds')

#initial conductivity choices for wadsleyite layer
pide_obj.set_mineral_conductivity_choice(rwd_wds = 5) #Dai and Karato 2009, Wadsleyite model
pide_obj.set_mineral_conductivity_choice(cpx = 8) #Xu
pide_obj.set_mineral_conductivity_choice(perov = 0)
pide_obj.set_mineral_conductivity_choice(garnet = 0)

pide_obj.list_transition_zone_water_partitions_solid('garnet')
pide_obj.list_transition_zone_water_partitions_solid('perov')
pide_obj.list_transition_zone_water_partitions_solid('cpx')

pide_obj.set_mantle_water_solubility(cpx = 2, perov = 0, garnet = 2)

for index in range(3,len(cond_arrays)):

	print('STARTED CALCULATIONS')

	if index == 6:
		
		pide_obj.set_mineral_conductivity_choice(rwd_wds = 1) #Dai and Karato 2009, Wadsleyite model
		
	pide_obj.set_temperature(T_global[:,idx_depth_list[index]] * np.ones(len(cond_arrays[index])))
	pide_obj.set_mantle_transition_zone_water_partitions(garnet = 0, perov = 0, cpx = 0)
	pide_obj.set_pressure(p[index])
	pide_obj.revalue_arrays()
	if index < 6:
		pide_obj.set_composition_solid_mineral(cpx = cpx[index], garnet = garnet[index], perov = pv[index], rwd_wds = wad[index])
	else:
		pide_obj.set_composition_solid_mineral(cpx = cpx[index], garnet = garnet[index], perov = pv[index], rwd_wds = ring[index])
	pide_obj.set_solid_phs_mix_method(method = 1) #hashin-shtrikman lower-bound
	tz_solubility = pide_obj.calculate_transition_zone_water_solubility(method = 'array') #calculating solubility at the given temperature
	
	c_list, residual_list = conductivity_solver_single_param(object = pide_obj, cond_list = cond_arrays[index], param_name = 'bulk_water', upper_limit_list = tz_solubility,
		lower_limit_list= np.zeros(len(tz_solubility)), search_start = 10, acceptence_threshold = 0.5, transition_zone = True, num_cpu = 6)	
	
	write_slice(p = pide_obj.p, t = pide_obj.T, water = c_list, cond = cond_arrays[index], lat = lat_arrays, lon = lon_arrays, depth = depth_p[index], filename = 'TZ_Water_Content_' + str(index) + '.csv')

	



	







