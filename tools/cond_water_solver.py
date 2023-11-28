import numpy as np
import matplotlib.pyplot as plt

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')

sys.path.append(core_path_ext)
import SEL
from inversion import conductivity_solver_single_param

temp = np.arange(600,1300,5) #setting up temperature array
index_list = range(0,len(temp))
sel_obj = SEL.SEL() #creating the initial object

sel_obj.set_composition_solid_mineral(ol = 0.65, opx = 0.2, cpx = 0.1, garnet = 0.05)
sel_obj.set_parameter('ti_ol', 0.1)
sel_obj.set_bulk_water(100)
sel_obj.calculate_bulk_mantle_water_solubility(method = 'array')



sys.exit()
conductivity_solver_single_param(object = sel_obj, temp_array = temp, p_array = p, param_name = 'bulk_water', upper_limit_list = 1000, lower_limit_list= 0 , increment =50, num_cpu = 1)
