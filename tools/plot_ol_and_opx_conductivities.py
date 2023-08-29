#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEEL')

sys.path.append(core_path_ext)

import SEEL
import numpy as np
import matplotlib.pyplot as plt

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

temp = np.arange(600,1300,5) #setting up temperature array
a = SEEL.SEEL() #creating the initial object
a.set_composition_solid_mineral(ol_frac = [1.0]) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)
a.set_param1_mineral(ol = np.ones(len(temp)) * 0)
a.set_mineral_water(ol = np.ones(len(temp)) * 100)
idx_ol = a.match_mineral_index('ol') #getting olivine index
list_olivine_models = a.list_mineral_econd_models('ol') #listing all olivine electrical conductivity methods

cond_olivine_lists = []

for i in range(0,len(list_olivine_models)):
	a.set_mineral_conductivity_choice(ol = i)
	cond = a.calculate_mineral_conductivity(method = 'array', min_idx= idx_ol)
	cond_olivine_lists.append(cond)
	
	
a.set_composition_solid_mineral(opx_frac = [1.0])	
a.set_mineral_water(opx = np.ones(len(temp)) * 100)
idx_opx = a.match_mineral_index('opx') #getting olivine index
list_opx_models = a.list_mineral_econd_models('opx') #listing all olivine electrical conductivity methods

cond_opx_lists = []

for i in range(0,len(list_opx_models)):
	a.set_mineral_conductivity_choice(opx = i)
	cond = a.calculate_mineral_conductivity(method = 'array', min_idx= idx_opx)
	cond_opx_lists.append(cond)


	
fig = plt.figure(figsize = (15,10))
ax = plt.subplot(121)
for i in range(0,len(list_olivine_models)):
	ax.plot(1e4/temp, cond_olivine_lists[i], label = list_olivine_models[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
ax.set_title('Olivine models with 100ppm water and dry conditions',fontsize = 9)
ax.set_ylim((1e-11,1))

ax2 = plt.subplot(122)
for i in range(0,len(list_opx_models)):
	ax2.plot(1e4/temp, cond_opx_lists[i], label = list_opx_models[i])
ax2.set_yscale('log')
ax2.set_xlabel('10000/T [$K^{-1}$]')
ax2.grid(which = 'both')
ax2.legend(fontsize = 6)
ax2.set_title('Orthopyroxene models with 100ppm water and dry conditions',fontsize = 9)
ax2.set_ylim((1e-11,1))

plt.show()



	

