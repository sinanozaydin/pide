#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEEL')

sys.path.append(core_path_ext)

#GENERAL mineral identifiers
#ol:olivine
#opx:orthopyroxene
#cpx:clinopyroxene
#garnet:garnet
#mica:mica
#amp:amphibole
#quartz:quartz
#plag:plagioclase
#kfelds:k-feldspar
#sulphide:sulphide
#graphite:graphite
#mixture:mixture
#other:other

#GENERAL rock identifiers
# granite:granite
# granulite:granulite
# sandstone:sandstone
# gneiss:gneiss
# amphibolite:amphibolite
# basalt:basalt
# mud:mudstone/shale
# gabbro:gabbro
# other_rock: other rocks


import SEEL
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(600,1300,5) #setting up temperature array
temp_melt = np.arange(1300,1800,5)
a = SEEL.SEEL() #creating the initial object

a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)



a.set_composition_solid_mineral(ol = 0.6,opx = 0.25,cpx = 0.1, garnet= 0.05) #setting composition
a.set_phase_interconnectivities(ol = 1, opx = 2, cpx = 4, gt = 4, melt = 1)
a.set_mineral_water(ol = 50, opx = 50*5.6, cpx = 100*5.6, gt = 50*0.8)


mixing_lists = a.list_phs_mix_methods()

cond_lists = []
for i in range(0,len(mixing_lists)):
	a.set_solid_phs_mix_method(method = i)
	cond = a.calculate_conductivity(method = 'array')
	cond_lists.append(cond)



lines = ['-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']

fig = plt.figure(figsize = (15,10))
ax = plt.subplot(121)
for i in range(0,len(mixing_lists)):
	ax.plot(1e4/temp, cond_lists[i], label = mixing_lists[i], linestyle = lines[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
ax.set_title('Wet Lherzolite with different solid phase mixing methods',fontsize = 9)
ax.set_ylim((1e-9,10))


a.set_temperature(temp_melt)
a.set_pressure(1.0)
a.set_composition_solid_mineral(ol = 0.6,opx = 0.25,cpx = 0.1, garnet= 0.05) #setting composition
a.set_phase_interconnectivities(ol = 1, opx = 2, cpx = 4, gt = 4, melt = 1)
a.set_mineral_water(ol = 50, opx = 50*5.6, cpx = 100*5.6, gt = 50*0.8)
a.set_solid_phs_mix_method(method = 0)
cond_dry_matrix_hT = a.calculate_conductivity(method = 'array')

melt_mixing_lists = a.list_phs_melt_fluid_mix_methods()

a.set_solid_phs_mix_method(method = 0)
a.set_melt_or_fluid_mode(mode = 'melt')
a.set_melt_properties(co2 = 10e4, water = 1e4)
melt_frac = 0.1
a.set_melt_fluid_frac(frac = melt_frac)


cond_melt_lists = []

for i in range(0,len(melt_mixing_lists)):
	
	
	a.set_solid_melt_fluid_mix_method(method = i)
	cond = a.calculate_conductivity(method = 'array')

	cond_melt_lists.append(cond)
	
cond_melt = a.calculate_melt_conductivity(method = 'array')

ax2 = plt.subplot(122)
ax2.plot(1e4/temp_melt, cond_dry_matrix_hT, label = 'Rock conductivity', linestyle = '-',linewidth = 4)
for i in range(0,len(mixing_lists)):
	ax2.plot(1e4/temp_melt, cond_melt_lists[i], label = melt_mixing_lists[i], linestyle = lines[i])
	
ax2.plot(1e4/temp_melt, cond_melt, label = 'Melt conductivity', linestyle = ':',linewidth = 4)

ax2.set_yscale('log')
ax2.set_xlabel('10000/T [$K^{-1}$]')
ax2.set_ylabel('Conductivity [S/m]')
ax2.grid(which = 'both')
ax2.legend(fontsize = 6)
ax2.set_title('Wet Lherzolite with ' +str(melt_frac*100) + '% melt fraction''different solid-melt mixing methods',fontsize = 9)
ax2.set_ylim((1e-5,100))
	

plt.show()