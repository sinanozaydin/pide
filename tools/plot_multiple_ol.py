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

temp = np.arange(600,2000,5) #setting up temperature array
a = pide.pide() #creating the initial object
a.set_composition_solid_mineral(ol = [1.0]) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)
a.set_param1_mineral(ol = np.ones(len(temp)) * 0)
a.set_mineral_water(ol = np.ones(len(temp)) * 20)
idx_ol = a.get_mineral_index('ol') #getting olivine index
list_olivine_models = a.list_mineral_econd_models('ol') #listing all olivine electrical conductivity methods

cond_olivine_lists = []

for i in range(0,3):
	if i == 0:
		a.set_mineral_conductivity_choice(ol = 5)
	elif i == 1:
		a.set_mineral_conductivity_choice(ol = 14)
	elif i == 2:
		a.set_mineral_conductivity_choice(ol = [5,14])
	cond = a.calculate_mineral_conductivity(method = 'array', min_idx= idx_ol)
	cond_olivine_lists.append(cond)
		
listies = ['-','-','--']
colories = ['r','b','k']
fig = plt.figure(figsize = (15,10))
ax = plt.subplot(111)
for i in range(0,len(cond_olivine_lists)):
	ax.plot(1e4/temp, cond_olivine_lists[i], label = str(i), linestyle = listies[i],color = colories[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
ax.set_title('Olivine models with 100ppm water and dry conditions',fontsize = 9)
ax.set_ylim((1e-11,1))


plt.show()



	

