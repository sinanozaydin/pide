import pide
import numpy as np
import unittest

temp = np.array([1000])
pres = np.array(2.0)

p_obj = pide.pide()

p_obj.set_temperature(temp)
p_obj.set_pressure(pres)

p_obj.set_composition_solid_mineral(ol = 0.6, opx = 0.3, cpx = 0.05,garnet = 0.05)
p_obj.set_bulk_water(200)

p_obj.mantle_water_distribute()

water_list = [p_obj.ol_water[0],p_obj.opx_water[0],p_obj.cpx_water[0],p_obj.garnet_water[0]]
print(water_list)

p_obj.set_melt_fluid_frac(0.02)
p_obj.set_bulk_water(2000.0)
p_obj.mantle_water_distribute()

water_list = [p_obj.ol_water[0],p_obj.opx_water[0],p_obj.cpx_water[0],p_obj.garnet_water[0], p_obj.melt_water[0]]
print(water_list)

p_obj.reset()

temp = np.array([2000])
pres = np.array(13.0)

p_obj.set_temperature(temp)
p_obj.set_pressure(pres)

p_obj.set_composition_solid_mineral(rwd_wds = 0.6, cpx = 0.2, garnet = 0.1, perov = 0.1)
p_obj.set_bulk_water(1000.0)
p_obj.transition_zone_water_distribute()

water_list = [p_obj.rwd_wds_water[0],p_obj.cpx_water[0],p_obj.garnet_water[0],p_obj.perov_water[0]]
print(water_list)