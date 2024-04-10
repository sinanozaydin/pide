#!/usr/bin/env python3

import os,sys

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
a.set_composition_solid_rock(granite = 1) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)
idx_granite = a.get_rock_index('granite') #getting granite index
list_granite_models = a.list_rock_econd_models('granite') #listing all granite conductivity models

cond_granite_lists = []

for i in range(0,len(list_granite_models)):
	a.set_rock_conductivity_choice(granite = i)
	cond = a.calculate_rock_conductivity(method = 'array', rock_idx= idx_granite)
	cond_granite_lists.append(cond)
	
	
a.set_composition_solid_rock(granulite = 1.0)	
idx_granulite = a.get_rock_index('granulite') #getting granulite index
list_granulite_models = a.list_rock_econd_models('granulite') #listing all granulite electrical conductivity methods

cond_granulite_lists = []

for i in range(0,len(list_granulite_models)):
	a.set_rock_conductivity_choice(granulite = i)
	cond = a.calculate_rock_conductivity(method = 'array', rock_idx= idx_granulite)
	cond_granulite_lists.append(cond)


	
fig = plt.figure(figsize = (15,10))
ax = plt.subplot(121)
for i in range(0,len(list_granite_models)):
	ax.plot(1e4/temp, cond_granite_lists[i], label = list_granite_models[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
ax.set_title('Granite models with dry conditions',fontsize = 9)
ax.set_ylim((1e-11,1))

ax2 = plt.subplot(122)
for i in range(0,len(list_granulite_models)):
	ax2.plot(1e4/temp, cond_granulite_lists[i], label = list_granulite_models[i])
ax2.set_yscale('log')
ax2.set_xlabel('10000/T [$K^{-1}$]')
ax2.grid(which = 'both')
ax2.legend(fontsize = 6)
ax2.set_title('Granulite models with dry conditions',fontsize = 9)
ax2.set_ylim((1e-11,1))

plt.show()



	

