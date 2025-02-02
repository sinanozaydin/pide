import os,sys
import time
import numpy as np
import pide
from pide.inversion import conductivity_metropolis_hastings_two_param, conductivity_solver_single_param
from pide.imaging.plot_distribution import plot_posterior_distribution_two_params, plot_posterior_distribution_heatmap_two_params
from pide.utils.utils import save_h5_files
import matplotlib.pyplot as plt
import seaborn as sns

# Observed value (target output)
T = np.array([1100.0])
P = np.array([2])

#Bounds for the search space.
water_min = np.zeros(len(T))
water_max = 2000 * np.ones(len(T))
phlg_min = 0.00 * np.ones(len(T))
phlg_max = 0.05 * np.ones(len(T))

p_obj = pide.pide()
p_obj.set_temperature(T)
p_obj.set_pressure(P)
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.38,mica = 0.02)
p_obj.set_param1_mineral(mica = 0.48)
p_obj.set_solid_phs_mix_method(2) #Upper




"""
p_obj.set_bulk_water(800)
p_obj.mantle_water_distribute()


cond = p_obj.calculate_conductivity()
print(cond)
import ipdb
ipdb.set_trace()
sys.exit()
"""

cond_external = [1e-2]
initial_water = 100
initial_phlg = 0.02
initial_params = [[initial_water, initial_phlg]]

sigma = 1e-2 * np.ones(len(T))#in log
n_iterations = 500000
proposal_std = [100,0.03]
burning = 10000

for i in range(1):
	start_time = time.time()
	samples, acceptance_rates, misfits, samples_all, misfits_all = conductivity_metropolis_hastings_two_param(object = p_obj, cond_list = cond_external, initial_params = initial_params,param_name_1 = 'bulk_water',
	param_name_2= "mica_frac", upper_limits = (water_max,phlg_max),
		lower_limits = (water_min,phlg_min), sigma_cond = sigma,proposal_stds=proposal_std
		,n_iter = n_iterations, burning = burning, transition_zone = False,num_cpu = 1)
	
	end_time = time.time()
	print(f'Time passed for processing: {end_time-start_time} seconds')
	
	water_samples = samples[0][:, 0]
	phlg_samples = samples[0][:, 1]

	water_samples_all = samples_all[0][:, 0]
	phlg_samples_all = samples_all[0][:, 1]
	
	print(acceptance_rates)
	
	plot_posterior_distribution_two_params(data_param_1 = water_samples,data_param_2 = phlg_samples, file_name = f"{i}_distr.png",save = True)
	plot_posterior_distribution_heatmap_two_params(data_param_1 = water_samples,
	data_param_2 = phlg_samples, param_1_min = 0, param_1_max = 2000,
	param_2_min = 0, param_2_max = 0.05,file_name = f"{i}_solution.png",save = True)
	
import ipdb
ipdb.set_trace()



