import numpy as np
import matplotlib.pyplot as plt

import pide
from pide.inversion import conductivity_solver_single_param
from pide.geodyn.geotherm import calculate_hasterok2011_geotherm

moho = 38 #km
max_depth = 250

temp, depth, pres, idx_LAB = calculate_hasterok2011_geotherm(SHF = 36, T_0 =25.0,max_depth = max_depth, moho = moho)

pide_obj = pide.pide() #creating the initial object
pide_obj.set_temperature(temp)

pide_obj.set_pressure(pres)
pide_obj.set_composition_solid_mineral(ol = 0.65, opx = 0.2, cpx = 0.1, garnet = 0.05)
# pide_obj.set_phase_interconnectivities(ol = 1, opx = 2, cpx = 4, gt = 4,mica = 1)
pide_obj.set_solid_phs_mix_method(2)

pide_obj.list_mantle_water_partitions_solid('opx')
pide_obj.list_mantle_water_partitions_solid('cpx')
pide_obj.list_mantle_water_partitions_solid('garnet')
pide_obj.set_mantle_water_partitions(opx_ol = 3, cpx_ol = 4, garnet_ol = 0)

pide_obj.revalue_arrays()

cond_list_to_invert = np.ones(len(temp)) * 1e4
cond_list_to_invert[75:100] = 1000
cond_list_to_invert[100:150] = 100
cond_list_to_invert[150:] = 50
cond_list_to_invert = 1.0 / cond_list_to_invert #converting to conductivity
import time
a = time.time()
c_list, residual_list = conductivity_solver_single_param(object = pide_obj, cond_list = cond_list_to_invert, param_name = 'mica_frac', upper_limit_list = np.ones(len(temp)) * 0.1,
lower_limit_list= np.zeros(len(temp)), search_start = 0.01, acceptence_threshold = 0.001, num_cpu = 1)
b = time.time()

print(b-a)

import matplotlib.pyplot as plt

cond_calced = pide_obj.calculate_conductivity()

fig = plt.figure()
ax = plt.subplot(121)
# ax.plot(cond_list,object.T,label = 'data')
# ax.plot(cond_calced,object.T, label = 'calced')

ax.plot(pide_obj.ol_frac,depth,color = 'g')
ax.plot(pide_obj.ol_frac + pide_obj.opx_frac,depth, color = 'b')
ax.plot(pide_obj.ol_frac + pide_obj.opx_frac + pide_obj.cpx_frac,depth, color = 'cyan')
ax.plot(pide_obj.ol_frac + pide_obj.opx_frac + pide_obj.cpx_frac + pide_obj.garnet_frac,depth, color = 'r')
ax.plot(pide_obj.ol_frac + pide_obj.opx_frac + pide_obj.cpx_frac + pide_obj.garnet_frac + pide_obj.mica_frac,depth, color = 'k')
ax.set_xlim(0,1)
ax.invert_yaxis()


ax2 = plt.subplot(122)
ax2.plot(1.0/cond_list_to_invert,depth)
ax2.plot(1.0/cond_calced,depth, linestyle = '--')
ax2.set_xscale('log')
ax2.invert_yaxis()
# plt.savefig('2.png',dpi = 300)
plt.show()