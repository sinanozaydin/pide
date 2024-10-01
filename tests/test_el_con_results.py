import pide
import numpy as np
import ipdb

temp = np.array([1000])
pres = np.array(2.0)

p_obj = pide.pide()

p_obj.set_temperature(temp)
p_obj.set_pressure(pres)
p_obj.set_param1_mineral(mica = 0.2, plag = 0.1)

mineral_list = p_obj.list_available_minerals()
all_cond_list = []
for i in range(0,len(mineral_list)):

	list_models = p_obj.list_mineral_econd_models(mineral_list[i])
	mineral_cond_list = []
	for j in range(0,len(list_models)):
		
		eval('p_obj.set_mineral_conductivity_choice(' + mineral_list[i] + '=' + str(j) + ')')
		eval(f"p_obj.set_mineral_water({mineral_list[i]}= 50)")
		cond = eval(f"p_obj.calculate_mineral_conductivity(min_idx = '{mineral_list[i]}')")
		mineral_cond_list.append(cond[0])
	all_cond_list.append(mineral_cond_list)
print(all_cond_list)
