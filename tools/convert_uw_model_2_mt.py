#!/usr/bin/env python3

import os,sys
import h5py

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')
core_path_ext_2 = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL/sel_src')
sys.path.append(core_path_ext)
sys.path.append(core_path_ext_2)

import SEL
from material import Material
from model import Model
from geodyn.read_uw_model import *
from geodyn.interpolate_fields import interpolate_2d_fields

#setting up source folder of the files
source_folder = os.path.join('.','example_data','uwconversion')

#setting up filename folders for the h5 files.
temp_fnm = os.path.join(source_folder,'temperature-501.h5')
pstrain_fnm = os.path.join(source_folder,'projPlasticStrain-501.h5')
material_fnm = os.path.join(source_folder,'projMaterialField-501.h5')
stress_fnm = os.path.join(source_folder,'projStressField-501.h5')
melt_fnm = os.path.join(source_folder,'projMeltField-501.h5')
strain_rate_fnm = os.path.join(source_folder,'strainRateField-501.h5')
pressure_fnm = os.path.join(source_folder,'pressureField-501.h5')
mesh_fnm = os.path.join(source_folder,'mesh.h5')

py_start_fnm = os.path.join(source_folder,'Col262.py')

#reading h5 files related to 2D UW model
temp_data, pressure_data, mesh_data, material_data, pstrain_data, stress_data, melt_data, strain_rate_data = read_h5_files(temp_h5 = temp_fnm, pressure_h5 = pressure_fnm, mesh_h5 = mesh_fnm,
material_h5= material_fnm, strain_h5 = pstrain_fnm, stress_h5 = stress_fnm, melt_h5 = melt_fnm, strain_rate_h5= strain_rate_fnm)

#reading startup py file to get material order and properties.

material_names = read_uw_material_names_from_py_input(py_start_fnm)

#reading 2d mesh params mesh itself, mesh_centers mesh, array in x direction, array in y direction, borders of the mesh[max_x, min_x, max_y, min_y]
mesh, mesh_center, x_mesh, y_mesh, borders_mesh = setup_2d_mesh(mesh_data)

#getting material_array
material_array, air_material_idx = setup_material(material_data, material_names)


#getting strain array
pstrain_array = setup_uw_data_array_PROJ(pstrain_data)
stress_array = setup_uw_data_array_PROJ(stress_data)
melt_array = setup_uw_data_array_PROJ(melt_data)
strain_rate_array = setup_uw_data_array_PROJ(strain_rate_data)

temp_array = setup_uw_data_array_PROJ(temp_data) #mesh
pressure_array = setup_uw_data_array_PROJ(pressure_data) / 1e9 #converting to gigapascal

# plot_2D_underworld_Field(x_array = mesh[0], y_array = mesh[1], Field = material_array,cblimit_up = len(material_names), cblimit_down = 0, log_bool=False, cb_name = 'tab20b')
# plot_2D_underworld_Field(x_array = mesh[0], y_array = mesh[1], Field = pstrain_array, cblimit_up = 1e2, cblimit_down = 1e-2, log_bool=True, cb_name = 'binary')
# plot_2D_underworld_Field(x_array = mesh[0], y_array = mesh[1], Field = temp_array, cblimit_up = 1500, cblimit_down = 600, log_bool=False, cb_name = 'coolwarm', label = 'regular_array',plot_save = True)
# plot_2D_underworld_Field(x_array = mesh_center[0], y_array = mesh_center[1], Field = pressure_array, cblimit_up = 10, cblimit_down = 0, log_bool=False, cb_name = 'coolwarm')


#Converting larger arrays into mesh_center locations.
temp_array = interpolate_2d_fields(mesh,temp_array,mesh_center)
pressure_array = interpolate_2d_fields(mesh_center,pressure_array,mesh_center)
melt_array = interpolate_2d_fields(mesh,melt_array,mesh_center)
pstrain_array = interpolate_2d_fields(mesh,pstrain_array,mesh_center)
material_array = interpolate_2d_fields(mesh,material_array,mesh_center,method = 'nearest') #using nearest for the material for them to stay integers
strain_rate_array = interpolate_2d_fields(mesh_center,strain_rate_array,mesh_center,method = 'nearest') #using nearest for the material for them to stay integers

plot_2D_underworld_Field(x_array = mesh_center[0], y_array = mesh_center[1], Field = material_array,cblimit_up = len(material_names), cblimit_down = 0, log_bool=False, cb_name = 'tab20b',plot_save = True, label = 'material.png')

sel_object = SEL.SEL()
# list_garnet_models = sel_object.list_mineral_econd_models('garnet')
# list_cpx_models = sel_object.list_mineral_econd_models('cpx')
# list_ol_models = sel_object.list_mineral_econd_models('ol')
# list_opx_models = sel_object.list_mineral_econd_models('opx')

list_gabbro_models = sel_object.list_rock_econd_models('gabbro')
list_granite_models = sel_object.list_rock_econd_models('granite')
#Eclogite object
#index 17 garnet is Liu2022 Wet

#DEFINING THE OBJECTS with SEL.material
Eclogite_Object = Material(name = 'Eclogite_Object', material_index = 3, calculation_type = 'mineral', composition = {'garnet':0.6,'cpx':0.4}, interconnectivities = {'garnet':1, 'cpx':1.5}, 
el_cond_selections = {'garnet':17,'cpx':0}, water = {'garnet':50,'cpx':0}, xfe = {'garnet':0.2, 'cpx':0.5}, phase_mixing_idx = 0)


Lithospheric_Mantle_Object = Material(name = 'Lithospheric_Mantle_Object', material_index = 4, calculation_type = 'mineral', composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5}, 
el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0}, water_distr = True, water = {'bulk':100}, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, phase_mixing_idx = 0)

Asthenospheric_Mantle_Object = Material(name = 'Asthenospheric_Mantle_Object', material_index = 5,calculation_type = 'mineral', composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5}, 
el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0}, water_distr = True, water = {'bulk':100}, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, phase_mixing_idx = 0)

Mantle_Wedge_Object = Material(name = 'Mantle_Wedge_Object', material_index = 6, calculation_type = 'mineral', composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5}, 
el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0}, water = {'bulk':200}, water_distr = True, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, phase_mixing_idx = 0)

Oceanic_Sediment = Material(name = 'Oceanic_Sediment',material_index = 7, calculation_type = 'value', resistivity_medium = 50.0)

Oceanic_Upper_Crust = Material(name = 'Oceanic_Upper_Crust',material_index = 8, calculation_type = 'rock', composition = {'gabbro':1.0})

Continental_Sediments_UP_Object = Material(name = 'Continental_Sediments_UP_Object',material_index = 9, calculation_type = 'value', resistivity_medium = 100.0)

Continental_Upper_Crust_UP_Object = Material(name = 'Continental_Upper_Crust_UP_Object',material_index = 10, calculation_type = 'rock', composition = {'granite':1.0},interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0)

Continental_Lower_Crust_UP_Object = Material(name = 'Continental_Lower_Crust_UP_Object',material_index = 11, calculation_type = 'rock', composition = {'granulite':1.0},interconnectivities = {'granulite':1},
el_cond_selections = {'granulite':0}, phase_mixing_idx = 0)

Decollement_LP = Material(name = 'Decollement_LP',material_index = 12, calculation_type = 'value', resistivity_medium = 100.0)

Continental_Sediments_LP_Object = Material(name = 'Continental_Sediments_LP_Object', material_index = 13,calculation_type = 'value', resistivity_medium = 100.0)

Continental_Upper_Crust_LP_Object = Material(name = 'Continental_Upper_Crust_LP_Object', material_index = 14,calculation_type = 'rock', composition = {'granite':1.0},interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0)

Continental_Lower_Crust_LP_Object = Material(name = 'Continental_Lower_Crust_LP_Object', material_index = 15,calculation_type = 'rock', composition = {'granulite':1.0},interconnectivities = {'granulite':1},
el_cond_selections = {'granulite':0}, phase_mixing_idx = 0)

# #NotSure 
Decollement_UP = Material(name = 'Decollement_UP',material_index = 16, calculation_type = 'value', resistivity_medium = 100.0)

Fault = Material(name = 'Fault', material_index = 17, calculation_type = 'value', resistivity_medium = 100.0)

#creating material_object_list:
material_object_list = [Eclogite_Object,Lithospheric_Mantle_Object,Asthenospheric_Mantle_Object,Mantle_Wedge_Object,Oceanic_Sediment,Oceanic_Upper_Crust, Continental_Sediments_UP_Object,
Continental_Upper_Crust_UP_Object,Continental_Lower_Crust_UP_Object,Decollement_LP,Continental_Sediments_LP_Object,Continental_Upper_Crust_LP_Object,
Continental_Lower_Crust_LP_Object,Decollement_UP,Fault]

material_skip_list = [None, 5, 50,None,None,None, None,None,None,None,None,None,None,None,None]
# material_object_list = [Eclogite_Object, Oceanic_Upper_Crust,Mantle_Wedge_Object, Decollement_UP]

#creating model_object
mt_model_object = Model(material_list = material_object_list, material_array = material_array, T = temp_array, P = pressure_array, melt = melt_array, p_strain = pstrain_array, strain_rate = strain_rate_array, material_node_skip_rate_list = material_skip_list)
backgr_cond = mt_model_object.calculate_conductivity(type = 'background')


plot_2D_underworld_Field(x_array = mesh_center[0], y_array = mesh_center[1], Field = backgr_cond, cblimit_up = 1e3, cblimit_down = 1e-8, log_bool=True, cb_name = 'Spectral_r',plot_save = True,label = 'cond.png')