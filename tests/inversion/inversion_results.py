import numpy as np
import matplotlib.pyplot as plt

import pide
from pide.inversion import conductivity_solver_single_param
from pide.geodyn.geotherm import calculate_hasterok2011_geotherm

moho = 38 #km
max_depth = 250

T, depth, p, idx_LAB = calculate_hasterok2011_geotherm(SHF = 36, T_0 =25.0,max_depth = max_depth, moho = moho)

print(list(T))

p_obj = pide.pide() #creating the initial object
p_obj.set_temperature(T)
p_obj.set_pressure(p)

p_obj.set_composition_solid_mineral(ol = 0.65, opx = 0.2, cpx = 0.1, garnet = 0.05)
p_obj.set_solid_phs_mix_method(2) #Hashin-Shtrikman

p_obj.set_mantle_water_solubility(ol = 4,opx = 3, cpx = 0, garnet = 0)
p_obj.set_mantle_water_partitions(opx_ol = 3, cpx_ol = 4, garnet_ol = 0)
p_obj.revalue_arrays()

p_obj.set_parameter('ti_ol', 0.01)
max_water = p_obj.calculate_bulk_mantle_water_solubility(method = 'array')

cond_list_to_invert = np.ones(len(T)) * 1e4 #first 75 km, just creating a array to change the rest.
cond_list_to_invert[75:100] = 1000 #from 75 to 100 km
cond_list_to_invert[100:150] = 100 #from 100 to 150 km
cond_list_to_invert[150:] = 50 #from 150 to 250 km.
cond_list_to_invert = 1.0 / cond_list_to_invert #converting to conductivity

c_list, residual_list = conductivity_solver_single_param(object = p_obj, cond_list = cond_list_to_invert,
param_name = 'bulk_water', upper_limit_list = max_water, lower_limit_list= np.zeros(len(max_water)),
search_start = 30, acceptence_threshold = 0.5, num_cpu = 5)

print(c_list)