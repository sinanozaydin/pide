import numpy as np
import matplotlib.pyplot as plt

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')

sys.path.append(core_path_ext)
import SEL
from inversion import conductivity_solver_single_param


temp = np.arange(600,1300,5) #setting up temperature array
pres = np.ones(len(temp))
index_list = range(0,len(temp))
sel_obj = SEL.SEL() #creating the initial object
sel_obj.set_temperature(temp)
sel_obj.set_pressure(1.0)
sel_obj.set_composition_solid_mineral(ol = 0.65, opx = 0.2, cpx = 0.1, garnet = 0.05)
sel_obj.set_phase_interconnectivities(ol = 1, opx = 2, cpx = 4, gt = 4)
sel_obj.set_parameter('ti_ol', 0.1)

sel_obj.list_mantle_water_solubilities('ol')
sel_obj.list_mantle_water_solubilities('cpx')
sel_obj.list_mantle_water_solubilities('opx')
sel_obj.list_mantle_water_solubilities('garnet')
sel_obj.set_mantle_water_solubility(ol = 4,opx = 3, cpx = 0, garnet = 0)

sel_obj.list_mantle_water_partitions_solid('opx')
sel_obj.list_mantle_water_partitions_solid('cpx')
sel_obj.list_mantle_water_partitions_solid('garnet')
sel_obj.set_mantle_water_partitions(opx_ol = 3, cpx_ol = 4, garnet_ol = 0)

max_water = sel_obj.calculate_bulk_mantle_water_solubility(method = 'array')

cond_list = np.ones(len(max_water[0])) * 1e-3 #1000 ohm meter

conductivity_solver_single_param(object = sel_obj, cond_list = cond_list, param_name = 'bulk_water', upper_limit_list = max_water[0],
lower_limit_list= np.zeros(len(max_water[0])), search_start = 50, num_cpu = 1)
