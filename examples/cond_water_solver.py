import numpy as np
import matplotlib.pyplot as plt

import pide
from pide.inversion import conductivity_solver_single_param


temp = np.arange(600,1300,5) #setting up temperature array
pres = np.ones(len(temp))
index_list = range(0,len(temp))
pide_obj = pide.pide() #creating the initial object
pide_obj.set_temperature(temp)
pide_obj.set_pressure(3)
pide_obj.set_composition_solid_mineral(ol = 0.65, opx = 0.2, cpx = 0.1, garnet = 0.05)
pide_obj.set_phase_interconnectivities(ol = 1, opx = 2, cpx = 4, gt = 4)
# pide_obj.set_melt_fluid_frac(0.01)
pide_obj.set_parameter('ti_ol', 0.1)

pide_obj.list_mantle_water_solubilities('ol')
pide_obj.list_mantle_water_solubilities('cpx')
pide_obj.list_mantle_water_solubilities('opx')
pide_obj.list_mantle_water_solubilities('garnet')
pide_obj.set_mantle_water_solubility(ol = 4,opx = 3, cpx = 0, garnet = 0)

pide_obj.list_mantle_water_partitions_solid('opx')
pide_obj.list_mantle_water_partitions_solid('cpx')
pide_obj.list_mantle_water_partitions_solid('garnet')
pide_obj.set_mantle_water_partitions(opx_ol = 3, cpx_ol = 4, garnet_ol = 0)

pide_obj.revalue_arrays()

max_water = pide_obj.calculate_bulk_mantle_water_solubility(method = 'array')

cond_list = np.ones(len(max_water)) * 1e-3 #1000 ohm meter

c_list, residual_list = conductivity_solver_single_param(object = pide_obj, cond_list = cond_list, param_name = 'bulk_water', upper_limit_list = max_water,
lower_limit_list= np.zeros(len(max_water)), search_start = 10, acceptence_threshold = 0.5, num_cpu = 5)

import matplotlib.pyplot as plt
	
pide_obj.set_bulk_water(c_list)
pide_obj.mantle_water_distribute(method = 'array')
cond_calced = pide_obj.calculate_conductivity(method = 'array')
fig = plt.figure()
ax = plt.subplot(121)
# ax.plot(cond_list,object.T,label = 'data')
# ax.plot(cond_calced,object.T, label = 'calced')

ax.plot(c_list,pide_obj.T)
ax.set_ylim(np.amax(pide_obj.T),np.amin(pide_obj.T))
ax.plot(max_water,pide_obj.T)
ax.legend()
# plt.savefig('2.png',dpi = 300)
plt.show()