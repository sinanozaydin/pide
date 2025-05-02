import pide
import numpy as np

h2o = np.linspace(0,3000,200)
temp = np.ones(len(h2o)) * 1500
pres = np.ones(len(h2o)) * 3

p_obj = pide.pide()

p_obj.set_temperature(temp)
p_obj.set_pressure(pres)
p_obj.set_melt_fluid_conductivity_choice(melt = 2)
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.4)
p_obj.set_melt_fluid_frac(0.05)

p_obj.set_bulk_water(h2o)
p_obj.mantle_water_distribute()
p_obj.calculate_density_fluid(method = 'array')

print(p_obj.dens_melt_fluid)
import matplotlib.pyplot as plt

plt.plot(p_obj.dens_melt_fluid,p_obj.h2o_melt)
plt.show()
import ipdb
ipdb.set_trace()
