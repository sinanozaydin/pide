import numpy as np
import matplotlib.pyplot as plt

import os,sys,csv
import netCDF4

def read_csv(filename,delim):

	#Simple function for reading csv files and give out filtered output for given delimiter (delim)

	file_obj = open(filename,'rt',encoding = "utf8") #Creating file object
	file_csv = csv.reader(file_obj,delimiter = delim) #Reading the file object with csv module, delimiter assigned to ','
	data = [] #Creating empty array to append data

	#Appending data from csb object
	for row in file_csv:
		data.append(row)

	#Filtering data for None elements read.
	for j in range(0,len(data)):
		data[j] = list(filter(None,data[j]))
	data = list(filter(None,data))

	return data
	
def write_slice(p,t,water,cond, lat, lon, depth, filename):

	lines = ['Latitude, Longitude,Depth[km],P[GPa],T[K],Water[ppm], cond[Sm]\n']
	
	for i in range(0,len(p)):
		
		if lon[i] > 180.0
			lines.append(','.join((str(lat[i]), str(lon[i] - 360.0), str(depth), str(p[i]),str(t[i]),str(water[i]),str(cond[i]) + '\n')))
		else:
			lines.append(','.join((str(lat[i]), str(lon[i]), str(depth), str(p[i]),str(t[i]),str(water[i]),str(cond[i]) + '\n')))
		
	filesave = open(filename,'w')
	filesave.writelines(lines)
	filesave.close()
	
	print('Slice file is written at : ' + filename)

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')

sys.path.append(core_path_ext)
import SEL
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
				lon_arrays.append(lon_model[k])
			
	cond_arrays.append(np.array(cond_slice))
	
		
file_write = "cond_trial.csv"
filesave = open(file_write,'w')
filesave.writelines(line_trial)
filesave.close
	
cond_arrays = np.array(cond_arrays)
lat_arrays = np.array(lat_arrays)
lon_arrays = np.array(lon_arrays)


sel_obj = SEL.SEL()
sel_obj.list_mineral_econd_models('garnet')
sel_obj.list_mineral_econd_models('cpx')
sel_obj.list_mineral_econd_models('perov')
sel_obj.list_mineral_econd_models('rwd_wds')

#initial conductivity choices for wadsleyite layer
sel_obj.set_mineral_conductivity_choice(rwd_wds = 5) #Dai and Karato 2009, Wadsleyite model
sel_obj.set_mineral_conductivity_choice(cpx = 8) #Xu
sel_obj.set_mineral_conductivity_choice(perov = 0)
sel_obj.set_mineral_conductivity_choice(garnet = 0)

sel_obj.list_transition_zone_water_partitions_solid('garnet')
sel_obj.list_transition_zone_water_partitions_solid('perov')
sel_obj.list_transition_zone_water_partitions_solid('cpx')

sel_obj.set_mantle_water_solubility(cpx = 2, perov = 0, garnet = 2)


for index in range(0,1):

	if index == 6:
		
		sel_obj.set_mineral_conductivity_choice(rwd_wds = 1) #Dai and Karato 2009, Wadsleyite model
		
	sel_obj.set_temperature(t_p[index] * np.ones(len(cond_arrays[index])))
	sel_obj.set_mantle_transition_zone_water_partitions(garnet = 0, perov = 0, cpx = 0)
	sel_obj.set_pressure(p[index])
	sel_obj.revalue_arrays()
	if index < 6:
		sel_obj.set_composition_solid_mineral(cpx = cpx[index], garnet = garnet[index], perov = pv[index], rwd_wds = wad[index])
	else:
		sel_obj.set_composition_solid_mineral(cpx = cpx[index], garnet = garnet[index], perov = pv[index], rwd_wds = ring[index])
	sel_obj.set_solid_phs_mix_method(method = 1) #hashin-shtrikman lower-bound
	tz_solubility = sel_obj.calculate_transition_zone_water_solubility(method = 'array') #calculating solubility at the given temperature
	
	c_list, residual_list = conductivity_solver_single_param(object = sel_obj, cond_list = cond_arrays[index], param_name = 'bulk_water', upper_limit_list = tz_solubility,
		lower_limit_list= np.zeros(len(tz_solubility)), search_start = 10, acceptence_threshold = 0.5, transition_zone = True, num_cpu = 6)	
	
	write_slice(p = sel_obj.p, t = sel_obj.T, water = c_list, cond = cond_arrays[index], lat = lat_arrays, lon = lon_arrays, depth = depth_p[index], filename = 'TZ_Water_Content_' + str(index) + '.csv')
# fig = plt.figure()
# ax = plt.subplot(111)
# ax.scatter(lon_arrays,lat_arrays, c = np.log10(cond_arrays[0]), marker = 's', linewidth = 0.2, edgecolor = 'k',cmap = 'Spectral_r')
# plt.show()
	



	







