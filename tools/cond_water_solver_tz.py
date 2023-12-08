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
index_tz_end = 32

filename_phase_file = "/home/sinan/Desktop/Research/Transition_Zone_Ben/Phase_proportion.csv"
data_phases = read_csv(filename = filename_phase_file, delim = ',')

depth_of_interest = depth_model[index_tz_start:index_tz_end]

depth_p = []

for i in range(1,len(data_phases)):

	if (float(data_phases[i][1]) > 410.0) and (float(data_phases[i][1]) <= 660.0):
		
		if any(element < 1.5 for element in abs(float(data_phases[i][1]) - depth_of_interest)):
			depth_p.append(float(data_phases[i][1]))




