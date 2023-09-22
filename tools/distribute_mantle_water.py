#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')

sys.path.append(core_path_ext)

import SEL
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(600,1300,5) 
a = SEL.SEL() #creating the initial object
a.set_temperature(temp)
a.set_composition_solid_mineral(ol = 0.6, opx = 0.25, cpx = 0.1, garnet = 0.05) #setting composition
a.set_phase_interconnectivities(ol = 1, opx = 2, cpx = 4, gt = 4, melt = 1)
a.set_pressure(1.0)

opx_part_list = a.list_mantle_water_partitions_solid(mineral_name = 'opx')
cpx_part_list = a.list_mantle_water_partitions_solid(mineral_name = 'cpx')
garnet_part_list = a.list_mantle_water_partitions_solid(mineral_name = 'garnet')

idx_ol = a.get_mineral_index('ol') #getting olivine index
idx_opx = a.get_mineral_index('opx') #getting olivine index
idx_cpx = a.get_mineral_index('cpx') #getting olivine index
idx_garnet = a.get_mineral_index('garnet') #getting olivine index

a.set_mantle_water_partitions(opx_ol = 3, cpx_ol = 6, garnet_ol = 0)
a.set_bulk_water(np.linspace(100,1000,len(temp)))


a.mantle_water_distribute(method = 'array') #can extract a.ol_water, a.opx_water from this object now

fig = plt.figure(figsize = (15,10))
ax = plt.subplot2grid(211)
ax.plot(a.ol_water,temp-273.15, label = 'ol')
ax.plot(a.opx_water,temp-273.15, label = 'opx')
ax.plot(a.cpx_water,temp-273.15, label = 'cpx')
ax.plot(a.garnet_water,temp-273.15, label = 'garnet')
ax.plot(a.bulk_water,temp-273.15,linestyle = '--', linewidth = 3, label = 'Bulk', color = 'k')
ax.legend(fontsize = 8)
ax.grid()
ax.set_ylabel('Temperature [$C^{\circ}$]')
ax.set_xlabel('Water content [$H_2O$ wt. ppm]')

cond = a.calculate_conductivity(method = 'array')

#calculating dry assemblage
a.set_bulk_water(0)
cond_dry = a.calculate_conductivity(method = 'array')

ax2 = plt.subplot(212)
ax2.plot(cond,temp-273.15,label = 'Wet')
ax2.plot(cond_dry,temp-273.15,label = 'Dry')
ax2.set_xscale('log')
ax2.set_xlabel('Conductivity [S/m]')
ax2.set_ylabel('Temperature [$C^{\circ}$]')
ax2.grid()
ax2.legend()



a.set_melt_fluid_frac(0.02)
a.set_mantle_water_partitions(opx_ol = 3, cpx_ol = 6, garnet_ol = 0)
a.set_bulk_water(np.linspace(100,1000,len(temp)))
a.mantle_water_distribute(method = 'array')
cond_with_melt = a.calculate_conductivity

ax3 = plt.subplot(221)
ax3.plot(a.ol_water,temp-273.15, label = 'ol')
ax3.plot(a.opx_water,temp-273.15, label = 'opx')
ax3.plot(a.cpx_water,temp-273.15, label = 'cpx')
ax3.plot(a.garnet_water,temp-273.15, label = 'garnet')
ax3.plot(a.bulk_water,temp-273.15,linestyle = '--', linewidth = 3, label = 'Bulk', color = 'k')
ax3.plot(a.solid_water,temp-273.15,linestyle = '--', linewidth = 3, label = 'Solid')
ax3.plot(a.melt_water,temp-273.15,linestyle = '--', linewidth = 3, label = 'Solid')



ax4 = plt.subplot(222)

plt.show()