#Importing the neccesarry python libraries.
import os
import pide
from pide.inversion import conductivity_solver_single_param, conductivity_metropolis_hastings_two_param
import numpy as np
import matplotlib.pyplot as plt

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

temperature = np.ones(100) * 1500
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

cond_external = np.ones(len(temperature)) * 1e-1

#Now setting up the conductivity_solver_single_param
#we are using the parameter melt_fluid_mass_frac to estimate instead volumetric parameter melt_frac.
#This allows us to calculate the thermal expansion and water partitioning among melt and solid minerla mixture
#with more accuracy.
c_list, residual_list = conductivity_solver_single_param(object = p_obj, cond_list = cond_external,
param_name = 'melt_fluid_mass_frac', upper_limit_list = np.ones(len(temperature)), lower_limit_list= np.zeros(len(temperature)),
search_start = 0.01, acceptence_threshold = 0.05, num_cpu = 5)


fig = plt.figure()
ax = plt.subplot(111)
ax.plot(bulk_water_array,c_list)
plt.show()

import ipdb
ipdb.set_trace()