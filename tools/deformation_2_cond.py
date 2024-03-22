#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide')
core_path_ext_2 = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide/pide_src')
sys.path.append(core_path_ext)
sys.path.append(core_path_ext_2)

import pide
from geodyn.deform_cond import plastic_strain_2_conductivity
import numpy as np
import matplotlib.pyplot as plt

strain = np.logspace(-2,2,100)
strains = [1e-2,1e2]

strains_log = np.log10(strains)

low_cond_err = 10
mid_cond_err = 10
high_cond_err = 10

temp = np.array([1000]) #setting up temperature array
a = pide.pide() #creating the initial object
a.set_composition_solid_mineral(ol = [1.0]) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)
idx_ol = a.get_mineral_index('ol') #getting olivine index
cond_low = a.calculate_mineral_conductivity(method = 'array', min_idx = idx_ol)
idx_sulphide = a.get_mineral_index('sulphide')
cond_high = a.calculate_mineral_conductivity(method = 'array', min_idx = idx_sulphide)

conds = [cond_low,cond_high]

conds_log = np.log10(conds)

cond_decay = 0.5
strain_decay = 0.1
threshold_strain = 10

cond_strain = plastic_strain_2_conductivity(strain = strain, low_cond= cond_low[0], high_cond = cond_high[0], low_strain = strains[0], high_strain = strains[1],function_method = 'exponential',
strain_decay_factor = strain_decay, conductivity_decay_factor = cond_decay,strain_percolation_threshold = None, low_cond_err = low_cond_err, mid_cond_err = mid_cond_err,
high_cond_err = high_cond_err)

# cond_strain_threshold = plastic_strain_2_conductivity(strain = strain, low_cond= cond_low[0], high_cond = cond_high[0], low_strain = strains[0], high_strain = strains[1],function_method = 'exponential',
# strain_decay_factor = strain_decay, conductivity_decay_factor = cond_decay,strain_percolation_threshold = None)

cond_strain_threshold = plastic_strain_2_conductivity(strain = strain, low_cond= cond_low[0], high_cond = cond_high[0], low_strain = strains[0], high_strain = strains[1],function_method = 'exponential',
strain_decay_factor = strain_decay, conductivity_decay_factor = cond_decay,strain_percolation_threshold = None)

#calculating the same thing plastic_strain_2_conductivity does for plotting purposes (plotting the middle point)
cond_decay = cond_decay*-1
mid_strains_log = ((strains_log[1] - strains_log[0]) / 2.0) + strains_log[0]
mid_conds_log = ((conds_log[1] - conds_log[0]) / 2.0) + conds_log[0]
mid_conds_log = ((mid_conds_log - conds_log[0]) * cond_decay) + mid_conds_log
mid_strains_log = ((mid_strains_log - strains_log[0]) * strain_decay) + mid_strains_log
strains = 10**np.array([strains_log[0],mid_strains_log,strains_log[1]])
conds = 10**np.array([conds_log[0],mid_conds_log,conds_log[1]])
cond_err = np.array([conds[0][0]* (low_cond_err / 1e2), conds[1][0] * (mid_cond_err / 1e2), conds[2][0] * (high_cond_err / 1e2)])

fig = plt.figure(figsize = (10,5))
ax = plt.subplot(121)
ax.plot(strain,cond_strain)
ax.errorbar(strains,conds.T[0],yerr = cond_err,fmt = 'o')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('Accumulated Plastic Strain')
ax.set_ylabel('Conductivity [S/m]')
ax.grid(which = 'both')
ax2 = plt.subplot(122)
ax2.plot(strain,cond_strain_threshold)
ax2.plot(strains,conds,'o')
ax2.set_xlabel('Accumulated Plastic Strain')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.grid(which = 'both')
ax2.set_title('Case with the treshold value at ' +str(threshold_strain))
plt.show()
# plt.savefig('base2.png',dpi = 300)