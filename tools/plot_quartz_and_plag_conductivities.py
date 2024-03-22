#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide')

sys.path.append(core_path_ext)

import pide
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
a = pide.pide() #creating the initial object
a.set_composition_solid_mineral(plag = [1.0]) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)
a.set_param1_mineral(plag = np.ones(len(temp)) * 0)
a.set_mineral_water(plag = np.ones(len(temp)) * 0)
idx_plag = a.get_mineral_index('plag') #getting olivine index
list_plag_models = a.list_mineral_econd_models('plag') #listing all olivine electrical conductivity methods

cond_plag_lists = []

for i in range(0,len(list_plag_models)):
	a.set_mineral_conductivity_choice(plag = i)
	cond = a.calculate_mineral_conductivity(method = 'array', min_idx = idx_plag)
	cond_plag_lists.append(cond)
	
	
a.set_composition_solid_mineral(quartz = [1.0])	
idx_quartz = a.get_mineral_index('quartz') #getting quartz index
list_quartz_models = a.list_mineral_econd_models('quartz') #listing all quartz electrical conductivity methods
cond_quartz_lists = []

for i in range(0,len(list_quartz_models)):
	a.set_mineral_conductivity_choice(quartz = i)
	cond = a.calculate_mineral_conductivity(method = 'array', min_idx = idx_quartz)
	cond_quartz_lists.append(cond)

fig = plt.figure(figsize = (15,10))
ax = plt.subplot(121)
for i in range(0,len(list_plag_models)):
	ax.plot(1e4/temp, cond_plag_lists[i], label = list_plag_models[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
ax.set_title('Plagioclase models in dry conditions',fontsize = 9)
ax.set_ylim((1e-11,1))

ax2 = plt.subplot(122)
for i in range(0,len(list_quartz_models)):
	ax2.plot(1e4/temp, cond_quartz_lists[i], label = list_quartz_models[i])
ax2.set_yscale('log')
ax2.set_xlabel('10000/T [$K^{-1}$]')
ax2.grid(which = 'both')
ax2.legend(fontsize = 6)
ax2.set_title('Quartz models in dry conditions',fontsize = 9)
ax2.set_ylim((1e-11,1))

plt.show()



	

