from pide.material import Material
from pide.model import Model
from pide.utils.utils import sort_through_external_list
#For building a geotherm.
from pide.geodyn.geotherm import calculate_hasterok2011_geotherm
import numpy as np
import sys, ipdb,time
start = time.time()
Eclogite_Object = Material(name = 'Eclogite_Object', material_index = 3, calculation_type = 'mineral', composition = {'garnet':0.6,'cpx':0.4}, interconnectivities = {'garnet':1, 'cpx':1.5}, 
el_cond_selections = {'garnet':17,'cpx':0}, water = {'garnet':50,'cpx':0}, xfe = {'garnet':0.2, 'cpx':0.5}, solid_phase_mixing_idx = 0, deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None},top = 0.0, bottom = 10.0)

destination_obj = Material()
destination_obj_2 = Material()

Eclogite_Object.copy_attributes(destination_obj)
destination_obj.top = Eclogite_Object.bottom
destination_obj.bottom = 30.0
destination_obj.name = 'whatever'

destination_obj.copy_attributes(destination_obj_2)
destination_obj_2.top = 30.0
destination_obj_2.bottom = 100

T, depth, p, idx_LAB = calculate_hasterok2011_geotherm(SHF = 40, T_0 =25.0,max_depth = 250,moho = 38.0)

layers = Model(material_list = [Eclogite_Object,destination_obj,destination_obj_2], T = T, P = p, depth = depth)

layers.calculate_geothermal_block(type = 'conductivity')
end = time.time()

print(end-start)