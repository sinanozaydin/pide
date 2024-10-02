import pide
import numpy as np
import unittest

temp = np.array([1000,1250])
pres = np.array(2.0)

p_obj = pide.pide()

p_obj.set_temperature(temp)
p_obj.set_pressure(pres)

mineral_names_sol = ['ol','opx','cpx','garnet','rwd_wds']
p_obj.set_parameter('ti_ol',0.01)

for mineral_name in mineral_names_sol:
	
	list_sols = eval(f"p_obj.list_mantle_water_solubilities(mineral_name='{mineral_name}')")
	
	for j in range(0,len(list_sols)):
		eval(f"p_obj.set_mantle_water_solubility({mineral_name} = {j})")
		max_water = eval(f"p_obj.calculate_mineral_water_solubility('{mineral_name}')")
		print(max_water)
