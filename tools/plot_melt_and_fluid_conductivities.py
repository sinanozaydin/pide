#!/usr/bin/env python3

import os,sys

import pide
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(1200,1800,5)
a = pide.pide() #creating the initial object
a.set_melt_fluid_frac(1.0)
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1)
a.set_melt_properties(co2 = np.ones(len(temp)) * 0, water = np.ones(len(temp)) * 1e4, na2o = np.ones(len(temp)) * 0.1)

list_melt_models = a.list_melt_econd_models() #listing all melt electrical conductivity methods
list_fluid_models = a.list_fluid_econd_models() #listing all fluid electrical conductivity methods

cond_melt_lists = []
cond_fluid_lists = []

for i in range(0,len(list_melt_models)):
	a.set_melt_fluid_conductivity_choice(melt = i)
	cond = a.calculate_melt_conductivity(method = 'array')
	cond_melt_lists.append(cond)

temp_fluids = np.arange(600,1300,5) #setting up temperature array
a.set_temperature(temp_fluids) #settin temperature array in K
a.set_pressure(1.0)
a.set_fluid_properties(salinity = np.ones(len(temp_fluids)) * 0.1)

for i in range(0,len(list_fluid_models)):
	a.set_melt_fluid_conductivity_choice(fluid = i)
	cond = a.calculate_fluids_conductivity(method = 'array')
	cond_fluid_lists.append(cond)
	
fig = plt.figure(figsize = (15,10))
ax = plt.subplot(121)
ax0 = ax.twiny()
for i in range(0,len(cond_melt_lists)):
	ax.plot(1e4/temp, cond_melt_lists[i], label = list_melt_models[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
ax.set_title('Melt models with 100ppm water and dry conditions',fontsize = 9)
ax.set_ylim((1e-2,1000))
ax0.set_xlabel('T [K]')
ax.set_xlim((1e4/1850, 1e4/1150))
ax0.set_xlim(1850,1150)

ax2 = plt.subplot(122)
ax20 = ax2.twiny()
for i in range(0,len(cond_fluid_lists)):
	ax2.plot(1e4/temp_fluids, cond_fluid_lists[i], label = list_fluid_models[i])
ax2.set_yscale('log')
ax2.set_xlabel('10000/T [$K^{-1}$]')
ax2.grid(which = 'both')
ax2.legend(fontsize = 6)
ax2.set_title('Fluid models with 100ppm water and dry conditions',fontsize = 9)
ax2.set_ylim((1e-2,1000))
ax20.set_xlabel('T [K]')
ax2.set_xlim((1e4/1350, 1e4/550))
ax20.set_xlim(1350,550)

plt.show()
#plt.savefig('1.png',dpi = 300)


	

