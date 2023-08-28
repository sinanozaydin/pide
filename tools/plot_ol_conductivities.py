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
a.set_composition_solid_mineral(mica_frac = [1.0]) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)
a.set_param1_mineral(mica = np.ones(len(temp)) * 0.52)
a.set_mineral_water(mica = np.ones(len(temp)) * 11)
idx_mineral = a.match_mineral_index('mica') #getting olivine index
list_models = a.list_mineral_econd_models('mica') #listing all olivine electrical conductivity methods

cond_lists = []

for i in range(0,len(list_models)):
	a.set_mineral_conductivity_choice(mica = i)
	cond = a.calculate_mineral_conductivity(method = 'array', min_idx= idx_mineral)
	cond_lists.append(cond)
	
fig = plt.figure(figsize = (7,5))
ax = plt.subplot(111)
for i in range(0,len(list_models)):
	ax.plot(1e4/temp, cond_lists[i], label = list_models[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
# ax.set_ylim((1e-11,1))
plt.show()
# plt.savefig('0.png',dpi = 300)


	

