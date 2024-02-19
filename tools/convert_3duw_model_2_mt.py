#!/usr/bin/env python3

import os,sys
import h5py
import numpy as np

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')
core_path_ext_2 = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL/sel_src')
sys.path.append(core_path_ext)
sys.path.append(core_path_ext_2)

import SEL
from material import Material
from model import Model
from geodyn.read_uw_model import *
from geodyn.interpolate_fields import interpolate_3d_fields
from geodyn.plot_models import *
from geodyn.material_process import *
from geodyn.write_uw_model import write_3d_field_h5
from mt.mt_model_conversion import convert_2DModel_2_MARE2DEM

source_folder = "/media/sinan/Geodyn_HDD/Geodyn/ProtPA300r"

#setting up filename folders for the h5 files.
temp_fnm = os.path.join(source_folder,'temperature-0.h5')
pstrain_fnm = os.path.join(source_folder,'projPlasticStrain-0.h5')
material_fnm = os.path.join(source_folder,'projMaterialField-0.h5')
stress_fnm = os.path.join(source_folder,'projStressField-0.h5')
melt_fnm = os.path.join(source_folder,'projMeltField-0.h5')
strain_rate_fnm = os.path.join(source_folder,'strainRateField-0.h5')
pressure_fnm = os.path.join(source_folder,'pressureField-0.h5')
mesh_fnm = os.path.join(source_folder,'mesh.h5')

py_start_fnm = os.path.join(source_folder,'PullApt300r.py')

temp_data = read_h5_file(temp_fnm)
pressure_data = read_h5_file(pressure_fnm)
mesh_data = read_h5_file(mesh_fnm)
material_data = read_h5_file(material_fnm)
pstrain_data = read_h5_file(pstrain_fnm)
stress_data = read_h5_file(stress_fnm)
melt_data = read_h5_file(melt_fnm)
strain_rate_data = read_h5_file(strain_rate_fnm)

#reading startup py file to get material order and properties.

material_names = read_uw_material_names_from_py_input(py_start_fnm)


material_array, air_material_idx = setup_material(material_data, material_names)


#reading 2d mesh params mesh itself, mesh_centers mesh, array in x direction, array in y direction, borders of the mesh[max_x, min_x, max_y, min_y]
mesh, mesh_center, x_mesh, y_mesh, z_mesh, x_mesh_centers, y_mesh_centers, z_mesh_centers, borders_mesh = setup_3d_mesh(mesh_data)



#getting strain array
pstrain_array = setup_uw_data_array_PROJ(pstrain_data)
stress_array = setup_uw_data_array_PROJ(stress_data)
melt_array = setup_uw_data_array_PROJ(melt_data)
strain_rate_array = setup_uw_data_array_PROJ(strain_rate_data)
temp_array = setup_uw_data_array_PROJ(temp_data) #mesh
pressure_array = setup_uw_data_array_PROJ(pressure_data) / 1e9 #converting to gigapascal

#Converting larger arrays into mesh_center locations.
temp_array = interpolate_3d_fields((x_mesh,y_mesh,z_mesh),temp_array,(x_mesh_centers, y_mesh_centers, z_mesh_centers))
melt_array = interpolate_3d_fields((x_mesh,y_mesh,z_mesh),melt_array,(x_mesh_centers, y_mesh_centers, z_mesh_centers))
pstrain_array = interpolate_3d_fields((x_mesh,y_mesh,z_mesh),pstrain_array,(x_mesh_centers, y_mesh_centers, z_mesh_centers))
material_array = interpolate_3d_fields((x_mesh,y_mesh,z_mesh),material_array,(x_mesh_centers, y_mesh_centers, z_mesh_centers))

#forming the main pide object
sel_object = SEL.SEL()

#listing conductivity models
list_fluid_models = sel_object.list_fluid_econd_models()
list_melt_models = sel_object.list_melt_econd_models()
list_granite_models = sel_object.list_rock_econd_models('granite')
list_granulite_models = sel_object.list_rock_econd_models('granulite')
list_plag_models = sel_object.list_mineral_econd_models('plag')
list_quartz_models = sel_object.list_mineral_econd_models('quartz')
list_garnet_models = sel_object.list_mineral_econd_models('garnet')
list_opx_models = sel_object.list_mineral_econd_models('opx')
list_amp_models = sel_object.list_mineral_econd_models('amp')


#DEFINING THE OBJECTS with SEL.material
UpperCrustObject = Material(name = 'UpperCrustObject',material_index = 3, calculation_type = 'rock', composition = {'granite':1.0},interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

UpperCrustObject_b = Material(name = 'UpperCrustObject',material_index = 3, calculation_type = 'rock', composition = {'granite':1.0},
interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0,
melt_fluid_incorporation_method = 'value', melt_or_fluid = 'fluid', melt_fluid_cond_selection = 0,melt_fluid_content = 0.1,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

UpperCrustObject1 = Material(name = 'UpperCrustObject_1',material_index = 4, calculation_type = 'rock', composition = {'granite':1.0},interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

UpperCrustObject1_b = Material(name = 'UpperCrustObject_1',material_index = 4, calculation_type = 'rock', composition = {'granite':1.0},
interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0,
melt_fluid_incorporation_method = 'value', melt_or_fluid = 'fluid', melt_fluid_cond_selection = 0,melt_fluid_content = 0.1,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

UpperCrustObject2 = Material(name = 'UpperCrustObject',material_index = 5, calculation_type = 'rock', composition = {'granite':1.0},interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

UpperCrustObject2_b = Material(name = 'UpperCrustObject',material_index = 5, calculation_type = 'rock', composition = {'granite':1.0},
interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0,
melt_fluid_incorporation_method = 'value', melt_or_fluid = 'fluid', melt_fluid_cond_selection = 0,melt_fluid_content = 0.1,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

UpperCrustObject3 = Material(name = 'UpperCrustObject',material_index = 6, calculation_type = 'rock', composition = {'granite':1.0},interconnectivities = {'granite':1},
el_cond_selections = {'granite':0}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

UpperCrustObject3_b = Material(name = 'UpperCrustObject',material_index = 6, calculation_type = 'rock', composition = {'granite':1.0},
interconnectivities = {'granite':1},
el_cond_selections = {'granite':0,'other_rock':3}, phase_mixing_idx = 0,
melt_fluid_incorporation_method = 'value', melt_or_fluid = 'fluid', melt_fluid_cond_selection = 0,melt_fluid_content = 0.1,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

#felsic granulite
LowerCrustObject1 = Material(name = 'Continental_Lower_Crust_UP_Object',material_index = 7, calculation_type = 'mineral', composition = {'plag':0.45, 'garnet': 0.25,
'opx':0.05, 'amp':0.1, 'quartz': 0.15},
el_cond_selections = {'plag':1, 'garnet': 0, 'opx':0, 'amp':0, 'quartz': 7}, phase_mixing_idx = 1,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

LowerCrustObject1_b = Material(name = 'Continental_Lower_Crust_UP_Object',material_index = 7, calculation_type = 'mineral', composition = {'plag':0.45, 'garnet': 0.25,
'opx':0.05, 'amp':0.1, 'quartz': 0.15},
el_cond_selections = {'plag':1, 'garnet': 0, 'opx':0, 'amp':0, 'quartz': 7}, phase_mixing_idx = 1,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

LowerCrustObject2 = Material(name = 'Continental_Lower_Crust_UP_Object',material_index = 8, calculation_type = 'rock', composition = {'granulite':1.0},interconnectivities = {'granulite':1},
el_cond_selections = {'granulite':0}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

LowerCrustObject2_b = Material(name = 'Continental_Lower_Crust_UP_Object',material_index = 8, calculation_type = 'rock', composition = {'granulite':0.95, 'other_rock':0.05},
interconnectivities = {'granulite':1, 'other_rock':1},
el_cond_selections = {'granulite':0,'other_rock':3}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

MantleLithosphereObject = Material(name = 'Lithospheric_Mantle_Object', material_index = 9, calculation_type = 'mineral', composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5}, 
el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0}, water_distr = True, water = {'bulk':100}, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, phase_mixing_idx = 0,deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

MantleLithosphereObject_b = Material(name = 'Lithospheric_Mantle_Object', material_index = 9, calculation_type = 'mineral', composition = {'ol':0.62,'opx':0.24,'garnet':0.045,'cpx':0.045,'sulphide':0.05},
interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5,'sulphide':1}, 
el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0,'sulphide':0}, water_distr = True, water = {'bulk':100}, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

MantleObject = Material(name = 'Asthenospheric_Mantle_Object', material_index = 10,calculation_type = 'mineral', composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5}, 
el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0}, water_distr = True, water = {'bulk':100}, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, phase_mixing_idx = 0,deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})

MantleObject_b = Material(name = 'Asthenospheric_Mantle_Object', material_index = 10,calculation_type = 'mineral', composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5}, 
el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0}, water_distr = True, water = {'bulk':300}, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, phase_mixing_idx = 0,
deformation_dict = {'function_method':'exponential',
'conductivity_decay_factor':0.2, 'strain_decay_factor':0.2,'strain_percolation_threshold':None})







