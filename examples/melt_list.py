#Importing the neccesarry python libraries.
import os
import pide
import numpy as np
import matplotlib.pyplot as plt

colors = [
    "#e6194b",  # Red
    "#3cb44b",  # Green
    "#0082c8",  # Blue
    "#f58231",  # Orange
    "#911eb4",  # Purple
    "#46f0f0",  # Cyan
    "#f032e6",  # Magenta
    "#fabebe",  # Pink
    "#008080",  # Teal
    "#e6beff",  # Lavender
    "#aa6e28",  # Brown
    "#800000",  # Maroon
    "#aaffc3",  # Mint
    "#808000",  # Olive
    "#000080",  # Navy
    "#808080",  # Gray
    "#ffffff",  # White
    "#000000",  # Black
    "#a9a9a9",  # Dark Gray
    "#ff69b4"   # Hot Pink
]


p_obj = pide.pide()
melt_econd = p_obj.list_melt_econd_models()

temperature = np.arange(1000,2000,20)
pressure = np.ones(len(temperature)) * 2.0

p_obj.set_temperature(temperature)
p_obj.set_pressure(pressure)

fig = plt.figure(figsize = (10,10))
ax = plt.subplot(111)
# ax2 = plt.subplot(122)
melt_cond_dry = []
melt_cond_wet = []
p_obj.set_melt_properties(water = 0)
p_obj.set_melt_fluid_frac(0.01)
p_obj.set_composition_solid_mineral(ol = 1)
p_obj.set_solid_melt_fluid_mix_method(1)
import time
for i in range(len(melt_econd)):
	
	p_obj.set_melt_fluid_conductivity_choice(melt = i)
	mc = p_obj.calculate_conductivity()
	melt_cond_dry.append(mc)
	print(melt_econd[i])
	print('MELT',p_obj.melt_fluid_cond)
	print('BULK',p_obj.bulk_cond)
	print('##############')

	
	ax.plot(1e4/temperature, melt_cond_dry[i],label = melt_econd[i] + '_Dry',color = colors[i])

p_obj.set_melt_properties(water = 10000)

for i in range(len(melt_econd)):
	
	p_obj.set_melt_fluid_conductivity_choice(melt = i)
	mc = p_obj.calculate_conductivity()
	melt_cond_wet.append(mc)
	try:
		ax.plot(1e4/temperature, melt_cond_wet[i],label = melt_econd[i] + '_Wet',color = colors[i],linestyle = '--')
	except:
		import ipdb
		ipdb.set_trace()

import ipdb
ipdb.set_trace()
	
ax.set_yscale('log')
ax.legend()
ax.set_ylim(1e-5,10)


# ax2.set_yscale('log')
# ax2.set_ylim(1e-5,10)
plt.show()