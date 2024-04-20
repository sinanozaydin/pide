from pide.material import Material
from pide.utils.utils import sort_through_external_list
import numpy as np
import sys, ipdb

Eclogite_Object = Material(name = 'Eclogite_Object', material_index = 3, calculation_type = 'mineral', composition = {'garnet':0.6,'cpx':0.4}, interconnectivities = {'garnet':1, 'cpx':1.5}, 
el_cond_selections = {'garnet':17,'cpx':0}, water = {'garnet':50,'cpx':0}, xfe = {'garnet':0.2, 'cpx':0.5}, solid_phase_mixing_idx = 0, deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None},top = 0.0, bottom = 10.0)

destination_obj = Material()

Eclogite_Object.copy_attributes(destination_obj)
destination_obj.top = Eclogite_Object.bottom
destination_obj.bottom = 30.0
destination_obj.name = 'whatever'

material_list = [Eclogite_Object,destination_obj]
idx_list = [material_list[0].top,material_list[1].top]

material_list_2 = sort_through_external_list(idx_list,material_list)
ipdb.set_trace()
sys.exit()


temp = np.arange(500,1500,5)
c = Eclogite_Object.calculate_conductivity(T = temp, P = 3 * np.ones(len(temp)), melt= np.zeros(len(temp)))

print(c)