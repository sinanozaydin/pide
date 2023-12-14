#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')

sys.path.append(core_path_ext)

import SEL
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

temp = np.arange(873,1273,5) #setting up temperature array
a = SEL.SEL() #creating the initial object
a.set_composition_solid_mineral(rwd_wds = [1.0]) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(15.0)
a.set_param1_mineral(rwd_wds = np.ones(len(temp)) * 0)
a.set_mineral_water(rwd_wds = np.ones(len(temp)) * 740)
idx_rwd_wds = a.get_mineral_index('rwd_wds') #getting olivine index
list_wds_models = a.list_mineral_econd_models('rwd_wds') #listing all olivine electrical conductivity methods



cond_rwd_wds_lists = []

for i in range(0,len(list_wds_models)):
	a.set_mineral_conductivity_choice(rwd_wds = i)
	cond = a.calculate_mineral_conductivity(method = 'array', min_idx= idx_rwd_wds)
	cond_rwd_wds_lists.append(cond)
	
	
fig = plt.figure(figsize = (15,10))
ax = plt.subplot(111)
for i in range(0,len(list_wds_models)):
	ax.plot(1e4/temp, cond_rwd_wds_lists[i], label = list_wds_models[i])
ax.set_yscale('log')
ax.set_xlabel('10000/T [$K^{-1}$]')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax.legend(fontsize = 6)
ax.set_title('Rwd/Wds models with 100ppm water and dry conditions',fontsize = 9)
ax.set_ylim((1e-4,100))


plt.show()



	

