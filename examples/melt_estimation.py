#Importing the neccesarry python libraries.
import os
import pide
from pide.inversion import conductivity_solver_single_param, conductivity_metropolis_hastings_two_param
import numpy as np
import matplotlib.pyplot as plt
from pide.imaging.plot_distribution import plot_posterior_distribution_two_params, plot_posterior_distribution_heatmap_two_params

"""
#creating a pide object
p_obj = pide.pide()

#setting up and environment at 1300 K and 3 GPa
temperature = np.array([1500])
pressure = np.array([3.0])

p_obj.set_temperature(temperature)
p_obj.set_pressure(pressure)

#setting up a simple lherzolite composition
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.3, cpx = 0.05, garnet = 0.05)
#setting the mineral conductivity choices
p_obj.set_mineral_conductivity_choice(ol = 4, opx = 0, garnet = 0, cpx = 6)
p_obj.set_solid_phs_mix_method(1) #HS lower bound

#Setting up Sifre2014 as the conductivity choice
p_obj.set_melt_fluid_conductivity_choice(melt = 0)

#Setting up melt CO2
p_obj.set_melt_properties(co2 = 100)

#Setting up bulk water.
#Try to avoid setting up melt_h2o from set_melt properties, because this will not distribute water among minerals.
p_obj.set_bulk_water(300)
p_obj.mantle_water_distribute()

#Setting up the melt-solid mixture relationship to Tubes model
p_obj.set_solid_melt_fluid_mix_method(1)

#Setting up the conductivity value to estimate melt for
cond_external = [1e-1] #or 10 ohm m

#Now setting up the conductivity_solver_single_param
#we are using the parameter melt_fluid_mass_frac to estimate instead volumetric parameter melt_frac.
#This allows us to calculate the thermal expansion and water partitioning among melt and solid minerla mixture
#with more accuracy.
c_list, residual_list = conductivity_solver_single_param(object = p_obj, cond_list = cond_external,
param_name = 'melt_fluid_mass_frac', upper_limit_list = np.ones(len(temperature)), lower_limit_list= np.zeros(len(temperature)),
search_start = 0.01, acceptence_threshold = 0.5, num_cpu = 1)
"""


p_obj = pide.pide()

temperature = np.ones(100) * 1373
pressure = np.ones(100) * 3
bulk_water_array = np.linspace(100,3000,100)

p_obj.set_temperature(temperature)
p_obj.set_pressure(pressure)

#setting up a simple lherzolite composition
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.3, cpx = 0.05, garnet = 0.05)
#setting the mineral conductivity choices
p_obj.set_mineral_conductivity_choice(ol = 4, opx = 0, garnet = 0, cpx = 6)
p_obj.set_solid_phs_mix_method(1) #HS lower bound

#Setting up Sifre2014 as the conductivity choice
p_obj.set_melt_fluid_conductivity_choice(melt = 0)

#Setting up melt CO2
p_obj.set_melt_properties(co2 = 100)

#Setting up bulk water.
#Try to avoid setting up melt_h2o from set_melt properties, because this will not distribute water among minerals.
p_obj.set_bulk_water(bulk_water_array)
p_obj.mantle_water_distribute()

#Setting up the melt-solid mixture relationship to Tubes model
p_obj.set_solid_melt_fluid_mix_method(1)
p_obj.set_mantle_water_partitions(opx_melt = 1,cpx_melt = 1)

cond_external = np.ones(len(temperature)) * (1.0/20)

#Now setting up the conductivity_solver_single_param
#we are using the parameter melt_fluid_mass_frac to estimate instead volumetric parameter melt_frac.
#This allows us to calculate the thermal expansion and water partitioning among melt and solid minerla mixture
#with more accuracy.
melt_frac_solution, residual_list = conductivity_solver_single_param(object = p_obj, cond_list = cond_external,
param_name = 'melt_fluid_mass_frac', upper_limit_list = np.ones(len(temperature)), lower_limit_list= np.zeros(len(temperature)),
search_start = 0.01, acceptence_threshold = 0.05, num_cpu = 1)

fig = plt.figure()
ax = plt.subplot(111)
ax.plot(np.array(melt_frac_solution) * 1e2, bulk_water_array)
ax.set_xlabel('Melt Fraction (%)')
ax.set_ylabel('Bulk Water Content (ppm)')
ax.grid(which = 'both')
# plt.show()
plt.savefig('melt_frac_estimate_single.png')
"""
p_obj = pide.pide()

temperature = [1500]
pressure = [3]

p_obj.set_temperature(temperature)
p_obj.set_pressure(pressure)

#setting up a simple lherzolite composition
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.3, cpx = 0.05, garnet = 0.05)
#setting the mineral conductivity choices
p_obj.set_mineral_conductivity_choice(ol = 4, opx = 0, garnet = 0, cpx = 6)
p_obj.set_solid_phs_mix_method(1) #HS lower bound

#Setting up Sifre2014 as the conductivity choice
p_obj.set_melt_fluid_conductivity_choice(melt = 0)

#Setting up melt CO2
p_obj.set_melt_properties(co2 = 100)

#Setting up the melt-solid mixture relationship to Tubes model
p_obj.set_solid_melt_fluid_mix_method(1)

p_obj.set_mantle_water_partitions(opx_melt = 1,cpx_melt = 1)

res_external = 20

#Bounds for the search space.
water_min = np.zeros(len(temperature))
water_max = 3000.0 * np.ones(len(temperature))
melt_min = 0.00 * np.ones(len(temperature))
melt_max = 0.5 * np.ones(len(temperature))

cond_external = [1.0/res_external]
initial_water = 1000
initial_melt = 0.3
initial_params = [[initial_water, initial_melt]]

sigma = 0.01 * np.ones(len(temperature))#in log
n_iterations = 200000
proposal_std = [200,0.25]
burning = 10000

samples, acceptance_rates, misfits, samples_all, misfits_all = conductivity_metropolis_hastings_two_param(object = p_obj, cond_list = cond_external,
initial_params = initial_params,param_name_1 = 'bulk_water',
param_name_2= "melt_fluid_mass_frac", upper_limits = (water_max,melt_max),
	lower_limits = (water_min,melt_min), sigma_cond = sigma,proposal_stds=proposal_std
	,n_iter = n_iterations, burning = burning, transition_zone = False,num_cpu = 1,adaptive_alg = True,
	step_size_limits = [25000,0.5])

water_samples = samples[0][:, 0]
melt_samples = samples[0][:, 1]

water_samples_all = samples_all[0][:, 0]
melt_samples_all = samples_all[0][:, 1]

plot_posterior_distribution_two_params(data_param_1 = melt_samples*1e2,data_param_2 = water_samples, file_name = f"melt_frac_distr.png",save = True)
plot_posterior_distribution_heatmap_two_params(data_param_1 = melt_samples*1e2,
data_param_2 = water_samples, param_1_min = 0, param_1_max = np.amax(melt_samples_all)*1e2,
param_2_min = 0, param_2_max = np.amax(water_samples_all),
param1_name = 'Melt Fraction (%)',
param2_name = "Water Content (ppm)",
file_name = f"melt_frac_solution.png",save = True)
"""

import ipdb
ipdb.set_trace()