from pide.material import Material
import numpy as np

Eclogite_Object = Material(name = 'Eclogite_Object', material_index = 3, calculation_type = 'mineral', composition = {'garnet':0.6,'cpx':0.4}, interconnectivities = {'garnet':1, 'cpx':1.5}, 
el_cond_selections = {'garnet':17,'cpx':0}, water = {'garnet':50,'cpx':0}, xfe = {'garnet':0.2, 'cpx':0.5}, solid_phase_mixing_idx = 0, deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

destination_obj = Material()

Eclogite_Object.copy_attributes(destination_obj)

temp = np.arange(500,1500,5)
c = Eclogite_Object.calculate_conductivity(T = temp, P = 3 * np.ones(len(temp)), melt= np.zeros(len(temp)))

print(c)