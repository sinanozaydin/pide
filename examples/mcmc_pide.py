import os,sys
import time
import numpy as np
import pide
from pide.inversion import conductivity_metropolis_hastings_two_param, conductivity_solver_single_param
from pide.imaging.plot_distribution import plot_posterior_distribution_two_params, plot_posterior_distribution_heatmap_two_params,plot_misfit_distribution_two_params
import matplotlib.pyplot as plt
import seaborn as sns

# Observed value (target output)
T = np.array([1500.0])
P = np.array([2])

#Bounds for the search space.
water_min = np.zeros(len(T))
water_max = 2000 * np.ones(len(T))
melt_min = 0.00 * np.ones(len(T))
melt_max = 0.05 * np.ones(len(T))

p_obj = pide.pide()
p_obj.set_temperature(T)
p_obj.set_pressure(P)
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.4)

"""
p_obj.set_bulk_water(34.21875)
p_obj.mantle_water_distribute()
cond = p_obj.calculate_conductivity()

import ipdb
ipdb.set_trace()
"""

cond_external = [1e-2]
initial_water = 100
initial_melt = 0.005
initial_params = [[initial_water, initial_melt]]

sigma = 1e-2 * np.ones(len(T))#in log
n_iterations = 300000
proposal_std = [100,0.03]
burning = 10000

c_list, residual_list = conductivity_solver_single_param(object = p_obj, cond_list = cond_external,
param_name = 'bulk_water', upper_limit_list = np.ones(len(T))* 1e4, lower_limit_list= np.zeros(len(T)),
search_start = 30, acceptence_threshold = 0.5, num_cpu = 1)

for i in range(1):
	start_time = time.time()
	samples, acceptance_rates, misfits, samples_all, misfits_all = conductivity_metropolis_hastings_two_param(object = p_obj, cond_list = cond_external, initial_params = initial_params,param_name_1 = 'bulk_water',
	param_name_2= "melt_fluid_mass_frac", upper_limits = (water_max,melt_max),
		lower_limits = (water_min,melt_min), sigma_cond = sigma,proposal_stds=proposal_std
		,n_iter = n_iterations, burning = burning, transition_zone = False,num_cpu = 1)
	end_time = time.time()
	print(f'Time passed for processing: {end_time-start_time} seconds')
	
	water_samples = samples[0][:, 0]
	melt_samples = samples[0][:, 1]

	water_samples_all = samples_all[0][:, 0]
	melt_samples_all = samples_all[0][:, 1]
	
	print(acceptance_rates)
	
	plot_posterior_distribution_two_params(data_param_1 = water_samples,data_param_2 = melt_samples, file_name = f"{i}_distr.png",save = True)
	plot_posterior_distribution_heatmap_two_params(data_param_1 = water_samples,
	data_param_2 = melt_samples, param_1_min = 0, param_1_max = 2000,
	param_2_min = 0, param_2_max = 0.05,file_name = f"{i}_solution.png",save = True)
	plot_misfit_distribution_two_params(data_param_1 = water_samples_all,
	data_param_2 = melt_samples_all, misfits= misfits_all[0], param_1_min = 0, param_1_max = 2000,
	param_2_min = 0, param_2_max = 0.05,file_name = f"{i}_misfit_map.png",save = False)
	
plt.plot(misfits_all[0])
plt.show()
plt.plot(acceptance_rates[0])
plt.show()
import ipdb
ipdb.set_trace()



